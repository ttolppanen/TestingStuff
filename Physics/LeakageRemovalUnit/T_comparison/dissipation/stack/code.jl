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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:gamma] = 27
    param[:t] = 400
    param[:d] = 3
    param[:L] = 5
    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")
    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
    param[:w] = 1000
    param[:U] = 50
    param[:J] = 1

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
    traj = 1000
    diss_rates = [5, 10, 15, 20]
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    state = param[:state].val

    msrop = param[:effect].val.msr
    feedback = param[:effect].val.fb

    op = singlesite_n(d, L, 1) + singlesite_n(d, L, 2)
    obs = (state, id) -> expval(state[id], op[id])

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)

    for diss_i in diss_rates
        out = []
        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = subspace_split(bosehubbard(d, L; w = dw, U = dU, J = param[:J]), ranges, perm_mat)
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs; effect!, out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "d = $diss_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1 + n2", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()