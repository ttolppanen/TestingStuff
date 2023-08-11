using Distributed
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

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    calc_order = zeroone_orderparameter(L)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measuresitesrandomly!(state, msrop, p)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    dst_y = r[1, :]
    # r = solvetrajectories(r_f, traj; paral = :distributed)
    # for p in prob
    #     effect!(state) = measuresitesrandomly!(state, msrop, p; kwargs...)
    #     r_f() = krylovevolve(state0, H, dt, t, k; effect!)
    #     r = solvetrajectories(r_f, traj; paral = :distributed)
    #     out = traj_mean(r)
    #     push!(ent_y, out[1, 1])
    #     push!(fluc_y, out[2, 1])
    # end
    return [LineInfo(0:dt:t, dst_y, traj, "L = $L")]
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 1,
        :J => 1,
        :d => 2,
        :measurement => "std_n_feedback_01",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :k => 6,
        :L => 12,
        :p => 0.06,
        :traj => 1000
    )
    dst_layers = krylov(param[:d], param[:L], param[:dt], param[:t], param[:traj], param[:p]; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
    println("Saving Data")
    path = joinpath(@__DIR__, "order_parameter_test", "L$(param[:L])_p$(param[:p])")
    # addtrajectories(path, ent_layers...)
    makeplot(path, dst_layers...; xlabel = "t", ylabel = "distance", param...)
end

@time f()
