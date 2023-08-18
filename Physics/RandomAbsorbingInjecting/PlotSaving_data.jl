using Distributed
using JLD2
using Statistics
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

@everywhere function alternating_measurement!(state::AbstractVector{<:Number}, msr_op, msr_prob::Real, abs_prob, abs_op, inj_op)
    for i in 1:length(msr_op)
        if rand(Float64) < msr_prob
            measuresite!(state, msr_op, i)
            if rand(Float64) < abs_prob
                state .= abs_op[i] * state
                normalize!(state)
            else
                state .= inj_op[i] * state
                normalize!(state)
            end
        end
    end
end

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

function krylov(d, L, dt, t, traj, msr_prob, prob; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    abs_op = [isodd(i) ? singlesite(n_bosons_projector(d, 0), L, i) : singlesite([1 0; 0 1], L, i) for i in 1:L]
    inj_op = [isodd(i) ? singlesite([1 0; 0 1], L, i) : singlesite(n_bosons_projector(d, 1), L, i) for i in 1:L]
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    # calc_order = zeroone_orderparameter(L)
    n = nall(d, L) ./ L
    calc_order(state) = expval(state, n)
    H = bosehubbard(d, L; w, J, U)

    dst_y = []
    for abs_prob in prob
        effect!(state) = alternating_measurement!(state, msrop, msr_prob, abs_prob, abs_op, inj_op)
        r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!, save_only_last = true)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        out = traj_mean(r)[1, 1]
        push!(dst_y, out)
    end
    return LineInfo(prob, dst_y, traj, "L = $L")
end

function f(p)
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 1,
        :J => 1,
        :d => 2,
        :measurement => "std_n_alternating_to_0101",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :k => 6,
        :traj => 200,
        :msr_prob => p
    )
    prob = 0.0:0.05:1.0
    lines = []
    for L in 8:2:12
        line = krylov(param[:d], L, param[:dt], param[:t], param[:traj], param[:msr_prob], prob; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(lines, line)
    end
    println("Saving Data")
    path = joinpath(@__DIR__, "absorbing_test_p$(p)")
    # addtrajectories(path, dst_layers...)
    makeplot(path, lines...; xlabel = "absorbing probability", ylabel = "<N>", param...)
end

@time f(0.02)
@time f(0.03)
@time f(0.04)
