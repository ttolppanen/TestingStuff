using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function measure_ends_effect!(msr_op::Vector, 
    dt, msr_rate, feedback::Vector)
    time = [rand() / msr_rate]
    out(state) = begin
        time[1] += dt

        if time[1] < 1 / msr_rate
            return state
        end
        time[1] -= 1 / msr_rate

        for site in [length(msr_op)]
            measuresite!(state, msr_op, site)
            state .= feedback[site] * state
            normalize!(state)
        end
        return state
    end
    return out
end

function initial_parameters(L)
    param = OrderedDict()
    param[:gamma] = 27
    param[:dt] = 0.02
    param[:t] = 400
    param[:d] = 3
    state = normalize(bosonstack(1, L, 1; d = param[:d]) + allzero(param[:d], L))
    param[:state] = ValueInfo(state, "|10..> + |00..>")
    param[:alg] = "krylov"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = ValueInfo(param[:w] .+ (rand(L) .- 0.5) * 10, "10")
    param[:U] = 50
    param[:dU] = ValueInfo(2 * param[:dw].val .- (2 * param[:w] - param[:U]), "10")
    param[:J] = 1

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(3)
    # -----------
    traj = 2000
    Ls = [3,4,5]
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        state = param[:state].val

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        ops = param[:state].val * param[:state].val'
        obs = (state) -> expval(state, ops)

        dt = param[:dt]; t = param[:t]; k = param[:k]
        pa_k = PA_krylov(state, k)
        pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)

        out = []
        r_f(out, pa_k) = begin
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            H = H_from_param(L, param)
            krylovevolve(state, H, dt, t, k, pa_k, obs; effect!, out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "initial state overlap", copy_code = false, param...)
    # addtrajectories(path, line)
end

function H_from_param(L, param)
    dw = param[:w] .+ (rand(L) .- 0.5) * 10
    dU = 2 * dw .- (2 * param[:w] - param[:U])
    H = bosehubbard(param[:d], L; w = dw, U = dU, J = param[:J])
    return H
end

@time f()