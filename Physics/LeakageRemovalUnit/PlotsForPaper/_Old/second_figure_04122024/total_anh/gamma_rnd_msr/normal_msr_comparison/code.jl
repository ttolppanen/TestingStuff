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
    param[:dt] = 0.01
    param[:t] = 200
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 10
    param[:U] = 50
    param[:J] = 1

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rnd msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    traj = 9000
    Ls = [2, 3]
    gammas = exp.(range(log(0.01), log(100), length = 101))
    # -----------
    lines = []
    for L in Ls
        println("L = $L")
        param = initial_parameters(L)
        d = param[:d]
        state = param[:state].val

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        anh = anhall(d, L)
        obs = (state) -> expval(state, anh)

        dt = param[:dt]; t = param[:t]
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, deepcopy(state))

        out = []
        for gamma in gammas
            r_f(out, work_vector) = begin
                dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
                dU = 2 * dw .- (2 * param[:w] - param[:U])
                H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
                U = exp(-1.0im * Matrix(H) * dt)
                effect! = measure_ends_effect!(msrop, dt, gamma, feedback)
                exactevolve(state, work_vector, U, dt, t, obs; effect!, out, save_only_last = true)
            end
            @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, param...)
    addtrajectories(path, lines...)
    xlog_plot(@__DIR__)
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