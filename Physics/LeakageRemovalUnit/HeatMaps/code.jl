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
    dt, msr_rate, feedback::Vector; skip_subspaces = [])
    time = [0.0]
    out(state, subspace_id) = begin
        time[1] += dt

        if time[1] < 1 / msr_rate
            return state, subspace_id::Int64
        end
        time[1] -= 1 / msr_rate

        if in(subspace_id, skip_subspaces)
            return state, subspace_id::Int64
        end
        current_subspace_id = subspace_id
        for site in [1, length(msr_op[current_subspace_id])]
            msr_outcome = measuresite!(state[current_subspace_id], msr_op[current_subspace_id], site)
            new_id, op = feedback[current_subspace_id][site][msr_outcome]
            if new_id == current_subspace_id
                state[new_id] .= op * state[current_subspace_id]
            else
                mul!(state[new_id], op, state[current_subspace_id])
            end
            normalize!(state[new_id])
            current_subspace_id = new_id
        end
        return state, current_subspace_id::Int64
    end
    return out
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:gamma] = 0.1
    param[:t] = 80
    param[:d] = 3
    param[:L] = 7
    param[:state] = ValueInfo(productstate(param[:d], [0,0,0,2,0,0,0]), "|0002000>")
    param[:alg] = "krylov BN-subspace"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = ValueInfo(param[:w] .+ (rand(param[:L]) .- 0.5) * 10, "10")
    param[:U] = 50
    param[:dU] = ValueInfo(2 * param[:dw].val .- (2 * param[:w] - param[:U]), "10")
    param[:J] = 1

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "msr_n ends -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()

    # -----------
    traj = 100
    # -----------

    d = param[:d]; L = param[:L];
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    n_op = [subspace_split(singlesite_n(d, L, i), ranges, perm_mat) for i in 1:L]
    obs = [(state, id) -> expval(state[id], n_op[i][id]) for i in 1:L]

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(length(obs), length(0:dt:t)), traj, pa_k)

    H = subspace_split(bosehubbard(d, L; w = param[:dw].val, U = param[:dU].val, J = param[:J]), ranges, perm_mat)
    r_f(out, pa_k) = begin 
        effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
        krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs...; effect!, out)
    end
    solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
    r = traj_mean(pa_out)
    heatmap(1:L, 0:dt:t, transpose(r), xlabel = "site", ylabel = "t", legend = :none)
end

function linearize_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    traj = collect(values(plotinfo.lines))[1].traj
    pl = plot(yaxis = :log10, 
    title = plotinfo.title * ", traj = $traj", 
    xlabel = "t*gamma", 
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    xlims = (0, 30),
    ylims = (10^-1, 1)
    )
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        gamma = p / plotinfo.parameters[:dt]
        plot!(pl, line.x .* gamma, line.y, label = line.tag)
    end
    path = joinpath(path, "linearized plot.png")
    savefig(pl, path)
    return pl
end

function lin_reg_slope(x, y) # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    return (A \ y)[1]
end

function combine_slope_plots()
    path = joinpath(@__DIR__, "L12")
    pl = plot(xlabel = "p", 
    ylabel = "Exp decay slopes", 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright
    )
    for U in [7, 9, 10, 11, 13]
        data_path = joinpath(path, "U$U")
        plotinfo = load(joinpath(data_path, "data.jld2"), "plotinfo")
        ps = []
        slopes = []
        for line in values(plotinfo.lines)
            p = parse(Float64, line.tag[5:end])
            gamma = p / plotinfo.parameters[:dt]
            x = line.x .* gamma
            y = log10.(line.y)
            line_start_index = Int(floor(length(x) / 4))
            line_end_index = Int(floor(length(x) * 3 / 4))
            slope = lin_reg_slope(x[line_start_index : line_end_index], y[line_start_index : line_end_index])
            push!(ps, p)
            push!(slopes, slope)
        end
        traj = collect(values(plotinfo.lines))[1].traj
        title!(plotinfo.title * ", traj = $traj")
        plot!(pl, ps, slopes, label = "U = $U")
    end
    path = joinpath(path, "combined slope plot test.png")
    savefig(pl, path)
    return pl
end

@time f()