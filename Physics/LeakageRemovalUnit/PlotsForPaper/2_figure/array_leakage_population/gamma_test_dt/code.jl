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
    param[:dt] = 0.1
    param[:t] = 200
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Ls = [3]
    gammas = exp.(x_ticks(log(30), log(200), 23))
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        anh = anhall(d, L)
        obs = (state) -> expval(state, anh)

        out = []
        for gamma in gammas
            dt = min(param[:dt], 1 / (4 * gamma))
            r_f() = begin
                dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
                dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
                H = bosehubbard(d, L; w = dw, U = dU, J = 1)
                Uop = exp(-1im * dt * Matrix(H))
                effect! = periodic_msr_feedback!_function(msrop, dt, gamma, feedback, [L])
                exactevolve(param[:state].val, Uop, dt, param[:t], obs; effect!, save_only_last = true)
            end
            @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, add_traj = true, param...)
    # addtrajectories(path, lines...)
    xlog_plot(@__DIR__)
end

function xlog_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
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
# xlog_plot(@__DIR__)
