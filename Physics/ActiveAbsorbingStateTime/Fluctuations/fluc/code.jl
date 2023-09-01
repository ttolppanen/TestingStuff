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

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    calc_order = zeroone_orderparameter(L)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measuresitesrandomly!(state, msrop, p)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
    r = solvetrajectories(r_f, traj; paral = :distributed)
    y = calc_fluctuations(r)

    return LineInfo(0:dt:t, y, traj, "p = $p")
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
        :L => 10,
        :k => 6
    )
    traj = 500
    prob = 0.01:0.02:0.19
    dst_lines = []
    for p in prob
        order_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, p; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "fluc")
    # addtrajectories(path, dst_lines...)
    makeplot(path, dst_lines...; xlabel = "t", ylabel = "Steady State Fluctuations", param...)
end

function fix_x_axis()
    plotinfo = load(joinpath(@__DIR__, "data.jld2"), "plotinfo")
    traj = [i.traj for i in values(plotinfo.lines)][1]
    pl = plot(title = plotinfo.title * ", traj = $traj", 
    xlabel = "t", 
    ylabel = plotinfo.ylabel,
    titlefontsize = 9,
    dpi = 300)
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y, label = line.tag)
    end
    path = joinpath(@__DIR__, "plot.png")
    savefig(pl, path)
end

# @time f()
fix_x_axis()
