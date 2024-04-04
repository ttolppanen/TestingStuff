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

function initial_parameters(L)
    param = OrderedDict()
    param[:dt] = 0.02
    param[:t] = 400
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(1, L, 1; d = param[:d]), "|100...>")
    param[:alg] = "krylov sbspc"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 20
    param[:U] = 50
    param[:dU] = 20
    param[:J] = 1

    return param
end

function f()
    param = initial_parameters(3)
    # -----------
    traj = 2000
    Ls = [3,4,5,6,7]
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        initial_id = find_subspace(param[:state].val, indices)
        state = subspace_split(param[:state].val, ranges, perm_mat)

        ops = subspace_split(param[:state].val * param[:state].val', ranges, perm_mat)
        obs = (state, id) -> expval(state[id], ops[id])

        dt = param[:dt]; t = param[:t]; k = param[:k]
        pa_k = PA_krylov_sub(state, k)
        pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)

        out = []
        r_f(out, pa_k) = begin
            H = H_from_param(L, param, ranges, perm_mat)
            krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs; out)
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

function H_from_param(L, param, ranges, perm_mat)
    dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
    dU = 2 * dw .- (2 * param[:w] - param[:U])
    H = bosehubbard(param[:d], L; w = dw, U = dU, J = param[:J])
    return subspace_split(H, ranges, perm_mat)
end

@time f()