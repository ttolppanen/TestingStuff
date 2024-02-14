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
    ps = [0.01, 0.02, 0.03, 0.04, 0.05]
    plots_per_p = 5
    U = 10
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

    plots = []
    added_bar = false
    for p in ps
        effect!(state, id) = random_measurement_feedback!(state, id, msrop, p, feedback; skip_subspaces = 1)
        for _ in 1:plots_per_p
            r = krylovevolve(state, initial_id, H, dt, t, k, obs...; effect!)
            if !added_bar
                added_bar = true
                push!(plots, heatmap(1:L, 0:dt:t, transpose(r)))
            else
                push!(plots, heatmap(1:L, 0:dt:t, transpose(r), xlabel = "site", ylabel = "t", title = "p = $p", legend = :none))
            end
        end
    end

    path = @__DIR__
    title = PlotAndSave.parameters_to_title(; U, param...)
    pl = plot(plots..., layout = (length(ps), plots_per_p), size = (1500, 1500), plot_title = title, plot_titlefont = 10)
    savefig(pl, joinpath(path, "heatmap.png"))
    # x = 1:L
    # y = 0:dt:t
    # # f(x, y) = r[x, Int(floor(y / 0.02)) + 1]
    # # pl = contour(x, y, f)
    # c = cgrad([:black, :yellow, :red, :purple, :white], [0.99, 1.0, 1.99, 2.0])
    # pl = heatmap(x, y, transpose(r), color=c)
    # plot(pl)
end

@time f()
#exit()