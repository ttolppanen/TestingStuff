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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:diss_rate] = 0.03
    param[:t] = 400
    param[:d] = 3
    param[:L] = 5
    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = param[:d]), "|100...>")
    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
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
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)

    op = subspace_split(singlesite_n(d, L, 1), ranges, perm_mat)
    obs = (state, id) -> expval(state[id], op[id])

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub_dd(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)
    
    a_relations = [-1 for _ in eachindex(ranges)]
    a_relations[1] = 0
    diss_op = [sqrt(param[:diss_rate]) * singlesite_a(d, L, L)]
    H_dd = diss_op[1]' * diss_op[1]
    diss_op = [subspace_split(op, ranges, perm_mat, a_relations) for op in diss_op]

    for dw_i in dws
        out = []
        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
            H .-= 1.0im / 2 .* H_dd
            H = subspace_split(H, ranges, perm_mat)
            krylovevolve(state, initial_id, H, diss_op, [], dt, t, k, pa_k, obs; out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()