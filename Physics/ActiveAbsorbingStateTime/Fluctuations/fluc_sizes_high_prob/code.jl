using Distributed
using JLD2
using Plots
using Statistics
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function calc_fluctuations(result)
    out = zeros(size(result[1][1, :]))
    st_value = mean([traj[1, end] for traj in result])
    for traj in result
        line = traj[1, :]
        out .+= (line .- st_value).^2
    end
    return out ./ length(result)
end

function zeroone_orderparameter(L)
    d = 2 # not well defined for higher d
    order_parameter = singlesite_n(d, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        if isodd(i)
            order_parameter .+= singlesite_n(d, L, i)
        else
            order_parameter .+= singlesite([1 0; 0 1], L, 1) .- singlesite_n(d, L, i)
        end
    end
    order_parameter ./= L
    return state -> expval(state, order_parameter)
end

function krylov(d, L, dt, t, traj, prob; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    calc_order = zeroone_orderparameter(L)
    H = bosehubbard(d, L; w, J, U)

    out = []
    num_of_steps = length(0:dt:t)
    for p in prob
        println("p = $p")
        effect!(state) = measuresitesrandomly!(state, msrop, p)
        r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        y = calc_fluctuations(r)
        push!(out, mean(y[Int(floor(num_of_steps / 2)) : num_of_steps]))
    end

    return LineInfo(prob, out, traj, "L = $L")
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 0.0,
        :J => 1,
        :d => 2,
        :measurement => "std_n_feedback_01",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :k => 6
    )
    traj = 500
    prob = 0.19:0.08:(0.19 * 4)
    sizes = [8, 10, 12]
    dst_lines = []
    for L in sizes
        order_line = krylov(param[:d], L, param[:dt], param[:t], traj, prob; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "fluc_sizes_high_prob")
    # addtrajectories(path, dst_lines...)
    makeplot(path, dst_lines...; xlabel = "p", ylabel = "Steady State Fluctuations", param...)
    plot_to_log(path)
    plot_to_log2(path)
end

function plot_to_log(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = "SSFluc * L",
    titlefontsize = 9, 
    dpi = 300,
    yaxis = :log10)
    for line in values(plotinfo.lines)
        L = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y * L, label = line.tag)
    end
    path = joinpath(path, "log_y_L.png")
    savefig(pl, path)
end

function plot_to_log2(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = "SSFluc * L",
    titlefontsize = 9, 
    dpi = 300,
    yaxis = :log10)
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y, label = line.tag)
    end
    path = joinpath(path, "log_y.png")
    savefig(pl, path)
end

@time f()
