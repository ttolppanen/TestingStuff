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
    dt, msr_rate, feedback::Vector)
    time = [rand() / msr_rate]
    out(state) = begin
        time[1] += dt

        if time[1] < 1 / msr_rate
            return state
        end
        time[1] -= 1 / msr_rate

        for site in [length(msr_op)]
            measuresite!(state, msr_op, site)
            state .= feedback[site] * state
            normalize!(state)
        end
        return state
    end
    return out
end

function initial_parameters(L)
    param = OrderedDict()
    param[:dt] = 100
    param[:t] = 5000
    param[:d] = 3
    param[:L] = 2
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:U] = ValueInfo(250, "250MHz")
    param[:w] = ValueInfo(7500 / param[:U].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:U].val, "100MHz")
    
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
    traj = 5000
    Js = [4, 5, 6] # In Mhz / 2pi
    gammas = range(0.0001, 0.001, length = 50)
    # -----------
    lines = []
    for J in Js
        d = param[:d]; L = param[:L]
        state = param[:state].val

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        anh = anhall(d, L)
        obs = (state) -> expval(state, anh)

        pa_work_vec = deepcopy(state)
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, pa_work_vec)

        out = []
        for gamma in gammas
            r_f(out, pa) = begin
                dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
                dU = 2 * dw .- (2 * param[:w].val - 1)
                H = bosehubbard(d, L; w = dw, U = dU, J = J / param[:U].val)
                U = exp(-1.0im * Matrix(H) * param[:dt])
                effect! = measure_ends_effect!(msrop, param[:dt], gamma, feedback)
                exactevolve(state, pa, U, param[:dt], param[:t], obs; effect!, out, save_only_last = true)
            end
            @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "J = $J"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, param...)
    # addtrajectories(path, lines...)
    # xlog_plot(@__DIR__)
end

function xlog_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.01, 0.03, 0.1, 5, 10, 20, 50], [0.01, 0.03, 0.1, 5, 10, 20, 50]),
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