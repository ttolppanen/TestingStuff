using Distributed
using JLD2
using Statistics
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function zeroone_orderparameter(L)
    d = 2 # not well defined for higher d
    order_parameter = singlesite_n(d, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        if isodd(i)
            order_parameter .+= singlesite_n(d, L, i)
        else
            order_parameter .+= singlesite([1 0; 0 1], L, 1) .- singlesite_n(d, L, i)
        end
    end
    order_parameter ./= L
    return state -> expval(state, order_parameter)
end

function return_times(result, dt, t; width = 0.01)
    out = Float64[]
    for traj in result
        did_not_return = true
        for i in 1:length(traj[1, :])
            time_step_value = traj[1, i]
            if abs(1 - time_step_value) <= width
                push!(out, i * dt)
                did_not_return = false
                break
            end
        end
        if did_not_return
            push!(out, t)
        end
    end
    return out
end

function krylov(d, L, dt, t, traj, prob; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    bra_zeroone = zeroone(d, L)'
    calc_dst(state) = abs(bra_zeroone * state)
    # calc_order = zeroone_orderparameter(L)
    H = bosehubbard(d, L; w, J, U)

    dst_y = []
    for p in prob
        effect!(state) = measuresitesrandomly!(state, msrop, p)
        r_f() = krylovevolve(state0, H, dt, t, k, calc_dst; effect!, save_only_last = true)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        out = mean(return_times(r, dt, t))
        push!(dst_y, out)
    end
    return LineInfo(prob, dst_y, traj, "L = $L")
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 120,
        :w => 0,
        :U => 1,
        :J => 1,
        :d => 2,
        :measurement => "std_n_feedback_01",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :k => 6,
        :traj => 200
    )
    prob = 0.01:0.005:0.08
    lines = []
    for L in 8:2:12
        line = krylov(param[:d], L, param[:dt], param[:t], param[:traj], prob; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(lines, line)
    end
    println("Saving Data")
    path = joinpath(@__DIR__, "return_mean_time")
    # addtrajectories(path, dst_layers...)
    makeplot(path, lines...; xlabel = "p", ylabel = "<T>", param...)
end

function combine_plots()
    plotinfo6 = load(joinpath(@__DIR__, "dst_prob_6", "data.jld2"), "plotinfo")
    plotinfo8 = load(joinpath(@__DIR__, "dst_prob", "data.jld2"), "plotinfo")
    line6 = plotinfo6.lines["L = 6"]
    line8 = plotinfo8.lines["L = 8"]
    path = joinpath(@__DIR__, "dst_prob_combined")
    makeplot(path, line6, line8; xlabel = "p", ylabel = "distance", plotinfo6.parameters...)
end

@time f()
# combine_plots()
