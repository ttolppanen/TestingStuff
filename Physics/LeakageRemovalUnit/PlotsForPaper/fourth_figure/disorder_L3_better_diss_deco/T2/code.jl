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
using SparseArrays

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
    param[:gamma] = 0.03
    param[:t] = 400
    param[:d] = 3
    param[:L] = 3

    state = complex(spzeros(param[:d]))
    state[1] = 1.0
    state[2] = 1.0
    normalize!(state)
    param[:state] = ValueInfo(state, "|0>+|1> + thermal state")
    param[:T] = 0.1
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34
    J = 5 * 10^6 * 2 * pi
    param[:J_value] = ValueInfo(J, "5MHz * 2pi")
    param[:T_scale] = ValueInfo(kb / (hbar * J), "kb/hbarJ")
    param[:diss] = 0.000106 # 300 mikrosecond
    param[:deco] = 0.000318 # 100 mikrosecond

    param[:alg] = "krylov sbspc"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1500
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
    dws = [5, 10, 15, 20]
    # -----------
    lines = []
    d = param[:d]; L = param[:L]

    msrop = param[:effect].val.msr
    feedback = param[:effect].val.fb

    sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
    ops = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, ops)

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_dd(allzero(d, L), k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)

    diss_op = [sqrt(param[:diss]) * singlesite_a(d, L, i) for i in 1:L]
    deco_op = [sqrt(param[:deco]) * singlesite_n(d, L, i) for i in 1:L]
    dd_op = vcat(diss_op, deco_op)
    H_dd = complex(spzeros(d^L, d^L))
    for op in dd_op
        H_dd .+= -1.0im / 2 * op' * op
    end
    
    for dw_i in dws
        out = []

        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J]) .+ H_dd

            t_T = param[:T]
            t_w = param[:w] * param[:J_value].val
            t_U = param[:U] * param[:J_value].val
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)

            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            krylovevolve(state, H, dd_op, dt, t, k, pa_k, obs; effect!, out)
        end
        r_f(zeros(1, length(0:dt:t)), pa_k) # this is necessary for some reason?
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "t", ylabel = "sigma_x1 env", copy_code = false, param...)
    addtrajectories(path, lines...)
end

@time f()