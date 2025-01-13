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
    param[:t] = 200
    param[:dt] = ValueInfo(0.1, "min(0.1, 1/10Γ)")
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy, with diss"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")
    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Ls = [3]
    gammas = exp.(x_ticks(log(100), log(200), 9))[2:end]
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        state = param[:state].val

        anh = anhall(d, L)
        obs = (state) -> expval(state, anh)

        t = param[:t]
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, PA_exact_dd(state))

        out = []
        for gamma in gammas
            diss_op = [sqrt(gamma) * singlesite_a(d, L, L)]
            H_dd = -1.0im / 2 * diss_op[1]' * diss_op[1]

            t_dt = min(param[:dt].val, 1 / (gamma * 10))
            r_f(out, pa) = begin
                dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
                dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
                H = bosehubbard(d, L; w = dw, U = dU, J = 1) .+ H_dd
                Uop = exp(-1im * t_dt * Matrix(H))
                exactevolve(state, Uop, diss_op, t_dt, t, pa, obs; out, save_only_last = true)
            end
            @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "diss", ylabel = "anh", copy_code = false, add_traj = true, param...)
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
