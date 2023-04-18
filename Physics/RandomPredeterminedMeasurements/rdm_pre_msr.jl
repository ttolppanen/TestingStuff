using Distributed
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution
using HDF5

function make_title(; kwargs...)
    out = ""
    for (key, value) in kwargs
        out *= "$key = $value, "
    end
    return out[1:end-2]
end

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end
function make_1010_msr_op(d, L)
    msrop = measurementoperators(nop(d), L)
    proj_op1 = n_bosons_projector(d, 1)
    proj_op2 = n_bosons_projector(d, 0)
    v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
    unitary_after_measurement!(v_p, msrop)
    return msrop
end

function make_projections_to_zeroes(d, L)
    out = []
    m = n_bosons_projector(d, 0)
    for i in 1:L
        push!(out, singlesite(m, L, i))
    end
    return out
end

function make_projections_to_zeroone(d, L)
    proj_op0 = n_bosons_projector(d, 0)
    proj_op1 = n_bosons_projector(d, 1)
    return [isodd(i) ? singlesite(proj_op0, L, i) : singlesite(proj_op1, L, i) for i in 1:L]
end

function make_n_corr_func(d, L)
    m = singlesite_n(d, L, 1) * singlesite_n(d, L, 2)
    for i in 2:L-1
        m .+= singlesite_n(d, L, i) * singlesite_n(d, L, i + 1)
    end
    m .= m ./ (L - 1)
    return state -> expval(state, m)
end

function projected_state_number(state)
    random_numb = rand(Float64)
    total_prob = 0
    for i in eachindex(state)
        total_prob += abs2(state[i])
        if random_numb < total_prob
            return i
        end
    end
end

function krylov(d, L, dt, t, k, traj, prob, proj_prob; kwargs...)
    state0 = zeroone(d, L)
    H = bosehubbard(d, L; kwargs...)
    msr_op = measurementoperators(nop(d), L)
    proj_op = make_projections_to_zeroone(d, L)
    out = zeros(length(prob), traj, L)
    # calc_ent(state) = entanglement(d, L, state, Int(floor(L/2)))
    # n = nall(d, L) ./ L
    # calc_n_tot(state) = expval(state, n)
    # calc_n_corr = make_n_corr_func(d, L)
    for i in eachindex(prob)
        effect!(state) = random_predetermined_measurement!(state, msr_op, prob[i], proj_op, proj_prob)
        r_f() = krylovevolve(state0, H, dt, t, k, projected_state_number; effect!, save_only_last = true)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        for r_i in eachindex(r)
            out[i, r_i, :] = digits(Int(r[r_i][1, 1]) - 1, base=2, pad = L) |> reverse
        end
    end
    h5open("feedback_rate_distributions.h5", "cw") do file
        L_key = "L=$L"
        if !haskey(file, L_key)
            create_group(file, L_key)
            attributes(file[L_key])["L"] = L
        end
        g = file[L_key]
        r_key = "r=$proj_prob"
        if !haskey(g, r_key)
            create_group(g, r_key)
            attributes(g[r_key])["feedback rate"] = proj_prob
        end
        g = g[r_key]
        g["data"] = out
        attributes(g["data"])["probabilities"] = collect(prob)
    end
    # res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent, calc_n_tot, calc_n_corr; paral = :distributed)
    # ent = res[1, :]
    # n_dens = res[2, :]
    # n_corr = res[3, :]
    # return [LineInfo(prob, ent, traj, "L = $L, r = $proj_prob"), LineInfo(prob, n_dens, traj, "L = $L, r = $proj_prob"), LineInfo(prob, n_corr, traj, "L = $L, r = $proj_prob")]
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 5,
        :J => 1,
        :d => 2,
        :measurement => "rdm_pdtmd_to_zeroone",
        :save_before_effect => false,
        :state => "zeroone",
        :alg => "krylov",
        :k => 6,
        :traj => 10
    )
    h5open("feedback_rate_distributions.h5", "w") do file
    end
    # ent_res = []
    # n_res = []
    # n_corr_res = []
    for r in [0.5, 1.0]
        for L in [6, 10]
            krylov(param[:d], L, param[:dt], param[:t], param[:k], param[:traj], 0.01:0.005:0.06, r; w = param[:w], U = param[:U], J = param[:J])
            # push!(ent_res, res[1])
            # push!(n_res, res[2])
            # push!(n_corr_res, res[3])
        end
    end
    # path = joinpath(@__DIR__, "out_01_r_higher_ent")
    # addtrajectories(path, ent_res...)
    # # makeplot(path, ent_res...; xlabel = "p", ylabel = "S", param...)
    # path = joinpath(@__DIR__, "out_01_r_higher_n")
    # addtrajectories(path, n_res...)
    # # makeplot(path, n_res...; xlabel = "p", ylabel = "N", param...)
    # path = joinpath(@__DIR__, "out_01_r_higher_n_corr")
    # addtrajectories(path, n_corr_res...)
    # # makeplot(path, n_corr_res...; xlabel = "p", ylabel = "N_corr", param...)
    println("DONE")
end

@time f()