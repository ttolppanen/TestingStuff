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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:diss_rate] = 0.03
    param[:t] = 400
    param[:d] = 3
    param[:L] = 5
    state = normalize(bosonstack(1, param[:L], 1; d = param[:d]) + allzero(param[:d], param[:L]))
    param[:state] = ValueInfo(state, "|10..> + |00..>")
    param[:alg] = "krylov sbspc"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:U] = 50
    param[:J] = 1
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
    state = param[:state].val

    sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
    ops = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, ops)

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_dd(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)
    
    diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    H_dd = diss_op[1]' * diss_op[1]

    for dw_i in dws
        out = []

        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
            H .-= 1.0im / 2 .* H_dd
            krylovevolve(state, H, diss_op, dt, t, k, pa_k, obs; out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "sigma_x1 env", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()