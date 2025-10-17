using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections
using HDF5

function save_data_to_hdf5(plt)
    path = joinpath(@__DIR__, "hdf5_data", "fig7_data.h5")
    file = h5open(path, "w")
    for s in plt.series_list
        label = String(s[:label])
        println(label)
        if label == "L = 3"
            label = "random feedback measurement"
        else
            continue
        end
        x = collect(s[:x])
        y = s[:y]
        sg = create_group(file, label)
        write(sg, "y", y)
        write(sg, "x", x)
    end
    close(file)
end

function add_attributes_to_hdf5_data()
    path = joinpath(@__DIR__, "hdf5_data", "fig7_data.h5")
    file = h5open(path, "r+")
    attrs(file)["DOI"] = "10.5281/zenodo.15516345"
    attrs(file)["description"] = "Figure 7 data."
    
    g_c = file["random feedback measurement"]

    attrs(g_c)["description"] = "Array leakage population as a function of random measurement rate. Rate is in the units of the coupling J. Random measurement here means that in each time-step of the simulation there is a probability of p = Gamma * dt to measure the last site."
    close(file)
end

function rdm_msr_feedback!(state::AbstractVector{<:Number}, msr_op::MsrOpMatrixType, msr_prob::Real, feedback::AbstractVector)
    site = 3
    if rand(Float64) < msr_prob
        measuresite!(state, msr_op, site)
        state .= feedback[site] * state
        normalize!(state)
    end
    return state
end

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters(L)
    param = OrderedDict()
    param[:dt] = ValueInfo(0.1, "min(0.05, 1/Γ)")
    param[:t] = 200
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rnd msr end -> fb -> |0>")
    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Ls = [3]
    gammas = exp.(x_ticks(log(10^-5), log(200), 200))
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
            dt = min(param[:dt].val, 1 / (8 * gamma))
            msr_prob = gamma * dt 
            r_f() = begin
                dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
                dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
                H = bosehubbard(d, L; w = dw, U = dU, J = 1)
                Uop = exp(-1im * dt * Matrix(H))
                effect!(state) = rdm_msr_feedback!(state, msrop, msr_prob, feedback)
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

# @time f()
# xlog_plot(@__DIR__)
save_data_to_hdf5(xlog_plot(@__DIR__))
add_attributes_to_hdf5_data()
