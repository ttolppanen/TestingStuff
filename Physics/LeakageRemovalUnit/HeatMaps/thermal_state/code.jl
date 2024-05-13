using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections
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
    param[:t] = 180
    param[:dt] = 0.02
    param[:gamma] = 0.03
    param[:d] = 3
    param[:L] = 5
    state = complex(spzeros(param[:d]))
    state[param[:d]] = 1.0
    param[:state] = ValueInfo(state, "|2> + thermal state")
    param[:T] = 1
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34
    J = 10 * 10^6
    param[:T_scale] = ValueInfo(kb / (hbar * J), "kb/hbarJ")
    param[:alg] = "krylov BN-subspace"
    param[:k] = 6
    param[:H] = "BH d_w fixed |2> energy"
    param[:w] = 1000
    param[:dw] = ValueInfo(param[:w] .+ (rand(param[:L]) .- 0.5) * 10, "10")
    param[:U] = 50
    param[:dU] = ValueInfo(2 * param[:dw].val .- (2 * param[:w] - param[:U]), "10")
    param[:J] = 1.5

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "msr_n ends -> fb -> 0")
    return param
end

function f()

    # -----------
    n_plots = 6
    traj = 500
    # -----------
    plots = []
    param = initial_parameters()
    for _ in 1:n_plots
        d = param[:d]; L = param[:L]

        H_diag = diag(bosehubbard(d, L - 1; w = param[:dw].val[2:end], U = param[:dU].val[2:end], J = 0))
        ts = thermal_state(H_diag, param[:T] * param[:T_scale].val)
        ts_sample = vec_from_density_op(ts)
        state = kron(param[:state].val, ts_sample)  

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        n_op = [singlesite_n(d, L, i) for i in 1:L]
        obs = [(state) -> expval(state, n_op[i]) for i in 1:L]

        dt = param[:dt]; t = param[:t]; k = param[:k]
        pa_k = PA_krylov(state, k)
        pa_args, pa_out = pre_alloc_threads(zeros(length(obs), length(0:dt:t)), traj, pa_k)

        H = bosehubbard(d, L; w = param[:dw].val, U = param[:dU].val, J = param[:J])

        r_f(out, pa_k) = begin
            ts_sample = vec_from_density_op(ts)
            state = kron(param[:state].val, ts_sample)
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            krylovevolve(state, H, dt, t, k, pa_k, obs...; effect!, out)
        end

        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        map = heatmap(1:L, 0:dt:t, transpose(r), 
        xlabel = "site", 
        ylabel = "t", 
        legend = :none,
        title = PlotAndSave.parameters_to_title(; L, param...),
        titlefont = 11)
        push!(plots, map)
    end

    # title = PlotAndSave.parameters_to_title(; param...)
    pl = plot(plots..., layout = (1, n_plots), size = (1500, 1500))
    savefig(pl, joinpath(@__DIR__, "heatmap.png"))
end

@time f()