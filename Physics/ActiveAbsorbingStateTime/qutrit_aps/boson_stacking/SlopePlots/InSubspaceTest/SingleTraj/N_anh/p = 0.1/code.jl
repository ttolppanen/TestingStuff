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
    param[:t] = 30
    param[:d] = 3
    param[:L] = 12
    param[:state] = ValueInfo(productstate(param[:d], [isodd(i) ? 2 : 0 for i in 1:param[:L]]), "|2020...>")
    param[:alg] = "krylov BN-subspace"
    param[:k] = 6
    param[:H] = "BH"
    param[:w] = 0
    param[:J] = 1

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 0; 0 0 0; 0 0 1]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "msr_n -> fb: 0 -> |0>, 1 -> |0>, 2 -> |2>")
    return param
end

function f()
    param = initial_parameters()

    # -----------
    p = 0.1
    U = 10
    traj = 50
    # -----------

    d = param[:d]; L = param[:L];
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)
    H = subspace_split(bosehubbard(d, L; w = param[:w], J = param[:J], U), ranges, perm_mat)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    n_op = [subspace_split(singlesite_n(d, L, i), ranges, perm_mat) for i in 1:L]
    obs = [(state, id) -> expval(state[id], n_op[i][id]) for i in 1:L]

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(length(obs), length(0:dt:t)), traj, pa_k)

    plots = []
    effect!(state, id) = random_measurement_feedback!(state, id, msrop, p, feedback; skip_subspaces = 1)
    r_f(out, pa_k) = krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs...; effect!, out)

    solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
    f_anh(n) = n .* (n .- 1)
    push!(plots, LineInfo(0:dt:t, vec(sum(traj_mean(pa_out) ./ L, dims = 1)), traj, "N/L"))
    push!(plots, LineInfo(0:dt:t, vec(sum(traj_mean(f_anh.(pa_out)) ./ L, dims = 1)), traj, "sum(N(N-1))/L"))

    path = joinpath(@__DIR__, "p = $p")
    makeplot(path, plots...; xlabel = "t", U, p, param...)
end

@time f()
#exit()