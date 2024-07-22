using Distributed
using JLD2
using Plots
using LinearAlgebra
using SparseArrays
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections
using DSP

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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:gamma] = 27
    param[:t] = 400
    param[:d] = 3
    param[:L] = 5
    state = normalize(bosonstack(1, param[:L], 1; d = param[:d]) + allzero(param[:d], param[:L]))
    param[:state] = ValueInfo(state, "|10..> + |00..>")
    param[:alg] = "krylov sbspc"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 10
    param[:U] = 50
    param[:J] = 1
    param[:J_val] = "5Mhz"

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    traj = 100
    diss_rates = [4, 8, 12, 16] ./ 5000
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    state = param[:state].val

    msrop = param[:effect].val.msr
    feedback = param[:effect].val.fb

    sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
    ops = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, ops)

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)
    
    diss_ops = [singlesite_a(d, L, i) for i in 1:L]

    for diss_i in diss_rates
        out = []
        diss_op_i = sqrt(diss_i) .* diss_ops
        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
            for op in diss_op_i
                H .-= 1.0im / 2 .* op' * op
            end

            msr_effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            effect!(state) = begin
                diss_deco_effect!(state, diss_op_i)
                normalize!(state)
                msr_effect!(state)
            end
            krylovevolve(state, H, dt, t, k, pa_k, obs; effect!, out, krylov_alg = :arnoldi)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "d = $diss_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "sigma_x1 env", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()