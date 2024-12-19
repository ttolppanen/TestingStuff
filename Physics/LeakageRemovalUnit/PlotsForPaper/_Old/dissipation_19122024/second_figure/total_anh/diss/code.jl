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
    param[:t] = 200
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "krylov with dissipation"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 10
    param[:U] = 50
    param[:J] = 1

    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    traj = 500
    Ls = [5]
    diss_rates = exp.(range(log(0.01), log(100), length = 101))
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        initial_id = find_subspace(param[:state].val, indices)
        ssplit(x) = subspace_split(x, ranges, perm_mat)
        state = ssplit(param[:state].val)

        anh = ssplit(anhall(d, L))
        obs = (state, id) -> expval(state[id], anh[id])

        dt = param[:dt]; t = param[:t]; k = param[:k]
        pa_k = PA_krylov_sub_dd(state, k)
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, pa_k)

        out = []
        a_relations = [-1 for _ in eachindex(ranges)]
        a_relations[1] = 0
        for gamma in diss_rates
            diss_op = [sqrt(gamma) * singlesite_a(d, L, L)]
            H_dd = diss_op[1]' * diss_op[1]
            diss_op = [subspace_split(op, ranges, perm_mat, a_relations) for op in diss_op]
            # for i in 2:L
            #     H_dd .+= diss_op[i]' * diss_op[i]
            # end
            r_f(out, pa_k) = begin
                dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
                dU = 2 * dw .- (2 * param[:w] - param[:U])
                H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
                H .-= 1.0im / 2 .* H_dd
                H = ssplit(H)
                krylovevolve(state, initial_id, H, diss_op, [], dt, t, k, pa_k, obs; out, save_only_last = true)
            end
            @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(diss_rates, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "diss", ylabel = "anh", copy_code = false, param...)
    addtrajectories(path, lines...)
    xlog_plot(@__DIR__)
end

function xlog_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.01, 0.03, 0.1, 5, 10, 20, 50, 100], [0.01, 0.03, 0.1, 5, 10, 20, 50, 100]),
    title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 7,
    legend = :outertopright
    # ylims = (10^-5, 1)
    )
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y, label = line.tag)
    end
    path = joinpath(path, "log plot.png")
    savefig(pl, path)
    return pl
end

@time f()