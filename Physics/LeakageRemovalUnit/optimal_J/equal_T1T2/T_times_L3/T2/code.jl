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

function initial_parameters(J)
    param = OrderedDict()
    param[:U] = ValueInfo(250, "250MHz")
    U = param[:U].val
    param[:gamma] = 3/2 * (J / U)^2 # this is 3/2 J^2/U, in units of U

    param[:t] = 1000 * U / 5 # J was 5Mhz
    param[:dt] = param[:t] / 10000
    param[:d] = 3
    param[:L] = 3

    state = complex(spzeros(param[:d]))
    state[1] = 1.0
    state[2] = 1.0
    normalize!(state)
    param[:state] = ValueInfo(state, "|0>+|1> + thermal state")

    param[:T] = 0.1
    param[:diss] = ValueInfo(0.005 / U, "5kHz")
    param[:deco] = ValueInfo(0.005 / U, "5kHz")
    
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
    param[:w] = ValueInfo(7500 / U, "7.5GHz")
    param[:dw] = ValueInfo(100 / U, "100MHz")

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]

    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(1)
    # -----------
    traj = 1000
    Js = 1:2:20
    # -----------
    lines = []
    for J_i in Js
        param = initial_parameters(J_i)
        d = param[:d]; L = param[:L]
        
        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb
        
        sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
        op = singlesite(sigma_x, L, 1)
        obs = (state) -> expval(state, op)

        dt = param[:dt]; t = param[:t];
        pa_e = QuantumTimeEvolution.PA_exact_dd(allzero(d, L))
        pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_e)

        diss_op = [sqrt(param[:diss].val) * singlesite_a(d, L, i) for i in 1:L]
        deco_op = [sqrt(param[:deco].val) * singlesite_n(d, L, i) for i in 1:L]
        H_dd = complex(spzeros(d^L, d^L))
        for op in vcat(diss_op, deco_op)
            H_dd .+= -1.0im / 2 * op' * op
        end

        out = []
        r_f(out, pa_e) = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - 1)
            H = bosehubbard(d, L; w = dw, U = dU, J = J_i / param[:U].val)

            t_T = param[:T]
            t_w = param[:w].val * param[:U].val * 10^6 * 2 * pi
            t_U = param[:U].val * 10^6 * 2 * pi
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)

            H .+= H_dd
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            U = exp(-1.0im * Matrix(H) * param[:dt])
            exactevolve(state, U, vcat(diss_op, deco_op), dt, t, pa_e, obs; effect!, out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "J = $J_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "sigma_x1 env", copy_code = false, param...)
    # addtrajectories(path, lines...)
end

@time f()