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

function krylov(d, L, dt, t, traj, msr_prob, abs_prob; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    abs_op = [singlesite(n_bosons_projector(d, 0), L, i) for i in 1:L]
    inj_op = [singlesite(n_bosons_projector(d, 1), L, i) for i in 1:L]
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    # calc_order = zeroone_orderparameter(L)
    n = nall(d, L) ./ L
    calc_order(state) = expval(state, n)
    H = bosehubbard(d, L; w, J, U)
        
    effect!(state) = alternating_measurement!(state, msrop, msr_prob, abs_prob, abs_op, inj_op)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!, save_only_last = false)
    r = solvetrajectories(r_f, traj; paral = :distributed)
    out = traj_mean(r)[1, :]

    return LineInfo(0:dt:t, out, traj, "msr_prob = $msr_prob")
end

function f(p)
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 1,
        :J => 1,
        :d => 2,
        :L => 10,
        :measurement => "std_n_alternating_to_0101",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :k => 6,
        :traj => 500,
        :abs_prob => p
    )
    prob = 0.01:0.04:0.21
    lines = []
    for msr_prob in prob
        line = krylov(param[:d], param[:L], param[:dt], param[:t], param[:traj], msr_prob, param[:abs_prob]; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(lines, line)
    end
    println("Saving Data")
    path = joinpath(@__DIR__, "time_no_structure_abs_prob$(param[:abs_prob])")
    # addtrajectories(path, dst_layers...)
    makeplot(path, lines...; xlabel = "measurement probability", ylabel = "<N/L>", param...)
end

@time f(0.1)
@time f(0.9)
