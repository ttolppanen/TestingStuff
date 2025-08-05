using JLD2
using Plots
using OrderedCollections
using PlotAndSave
using QuantumTimes
using LRUFigureColors
using LaTeXStrings


function TS_fit()
    path = joinpath(dirname(@__DIR__), "Ts")
    T_fit(path, 0.01)
end
function T1_fit()
    path = joinpath(dirname(@__DIR__), "T1")
    T_fit(path, 0.1)
end
function T2_fit()
    path = joinpath(dirname(@__DIR__), "T2")
    T_fit(path, 0.05, 0.7)
end

function ts_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "Ts", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "Array Leakage Population"
    )
    fit_lines = TS_fit()
    for line in values(plotinfo.lines)
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    
    savefig(pl, joinpath(path, "Ts.png"))
    savefig(pl, joinpath(path, "Ts.svg"))
    return pl
end

function t1_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T1", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "First Site Population",
    legend = :bottomleft
    )
    fit_lines = T1_fit()
    for line in values(plotinfo.lines)
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    savefig(pl, joinpath(path, "T1.png"))
    savefig(pl, joinpath(path, "T1.svg"))
    return pl
end

function t2_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T2", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = L"$\langle \sigma_x \rangle \textrm{envelope}$",
    legend = :none
    )
    fit_lines = T2_fit()
    for line in values(plotinfo.lines)
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end

    # normal decoherence
    # t = collect(values(plotinfo.lines))[1].x
    # deph_rate = parse(Float64, plotinfo.parameters[:deco][1:end - 3])
    # U = parse(Float64, plotinfo.parameters[:U][1:end - 3])
    # deph_rate = deph_rate / U
    # plot!(pl, t, x -> exp(-x * deph_rate * 0.5), label = "decoherence")

    line = collect(values(plotinfo.lines))[1]
    vline!(pl, [line.x[end] * 0.9])
    vline!(pl, [line.x[end] * 0.05])
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    savefig(pl, joinpath(path, "T2.png"))
    savefig(pl, joinpath(path, "T2.svg"))
    return pl
end



function make_fit_line(dir, cut_size)
    path = joinpath(dirname(@__DIR__), dir)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")

    lines = values(plotinfo.lines)
    t_times = []
    t_time_x = []
    for line in lines
        t_x = parse(Int64, line.tag[5:end])
        fit_line_start = floor(Int, length(line.x) * cut_size)
        t_line, t_time = T1(line.x[fit_line_start:end], line.y[fit_line_start:end])

        push!(t_times, t_time)
        push!(t_time_x, t_x)
    end
    return t_time_x, t_times
end

function combine_fit_plot()
    lines = TS_fit()
    t_time_xs = [parse(Float64, key[5:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    lines = T1_fit()
    t_time_x1 = [parse(Float64, key[5:end]) for key in keys(lines)]
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    t_time_x2 = [parse(Float64, key[5:end]) for key in keys(lines)]
    T2 = [val[2] for val in values(lines)]

    mega_cons = (10^6 * 2 * pi)
    data_path = joinpath(dirname(@__DIR__), "T1")
    J = get_param(data_path, :J) * 2 * pi

    pl = plot()


    plot!(pl, t_time_x1 * J / (2 * pi), T1 / J, label = L"$T_1$", color = color_t1)

    t1 = get_param(data_path, :diss)
    label1 = latexstring("T_1^{\\textrm{q}} = $(round(t1, digits=1)) \\ \\mu \\textrm{s}")
    hline!(pl, [t1], ls = :dash, label = "", color = color_t1)
    annotate!(pl, ((0.84, 0.91), text(label1, 7, "Computer Modern")))

    plot!(pl, t_time_x2 * J / (2 * pi), T2 / J, label = L"$T_2$", color = color_t2)

    pure_t2 = get_param(data_path, :deco)
    t2 = 1 / (1/pure_t2 + 1/(2*t1))
    label2 = latexstring("T_2^{\\textrm{q}} \\approx $(round(t2, digits=1)) \\ \\mu \\textrm{s}")
    hline!(pl, [t2], ls = :dash, label = "", color = color_t2)
    annotate!(pl, ((0.16, 0.40), text(label2, 7, "Computer Modern")))

    plot!(t_time_xs * J / (2 * pi), Ts / J,
    xlabel = L"$W / 2 \pi \ (\textrm{MHz})$",
    ylabel = L"$\textrm{Decay\ Times}\ (\mu \textrm{s})$",
    label = L"$T_\star$",
    legend = :topleft,
    color = color_ts,
    dpi = 300,
    titlefontsize = 9,
    ylims = [0, 17],
    xlims = [5, 100],
    yticks = 0:4:16,
    xminorticks = 4,
    yminorticks = 4,
    minorgrid = true
    )
    data_path = joinpath(dirname(@__DIR__), "T1")
    J = get_param(data_path, :J) * 2 * pi
    U = get_param(data_path, :U) * 2 * pi
    L = get_param(data_path, :L)
    t_y = t_small_scale(J, U, L)
    hline!(pl, [t_y], label = L"T_{\textrm{prop}}", ls = :dash, color = color_ts_short)

    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))

    lines = []
    traj = get_traj(data_path)
    push!(lines, LineInfo(t_time_xs, Ts, traj, "Ts"))
    push!(lines, LineInfo(t_time_x1, T1, traj, "T1"))
    push!(lines, LineInfo(t_time_x2, T2, traj, "T2"))
    plotinfo = loaddata(data_path)
    savedata(@__DIR__, lines...; xlabel = "Disorder", ylabel = "Survival Times", copy_code = false, plotinfo.parameters...)
    return pl
end

function ts_proportional_check()
    lines = TS_fit()
    t_time_xs = [parse(Float64, key[5:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    pl = plot(t_time_xs, Ts,
    xlabel = "J",
    ylabel = "Survival Times",
    label = "Ts",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9
    )
    L = 3
    U = 250
    f_ts(J) = 5/3*(L - 1) * 1 / (2 * (J/U)^2)
    plot!(pl, t_time_xs, f_ts, label = "5/3(L-1)U/(2J^2)")
    return pl
end

function print_last_val()
    lines = TS_fit()
    Ts = [val[2] for val in values(lines)]
    
    lines = T1_fit()
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    T2 = [val[2] for val in values(lines)]

    data_path = joinpath(dirname(@__DIR__), "T1")
    J = get_param(data_path, :J) * 2 * pi

    println("Ts = $(round(Ts[end] / J, digits = 2))")
    println("T1 = $(round(T1[end] / J, digits = 2))")
    println("T2 = $(round(T2[end] / J, digits = 2))")
end

#Plots.scalefontsizes(1.5)
# ts_plot()
# t1_plot()
# t2_plot()
combine_fit_plot()
# println("Trajectories: $(get_traj(joinpath(dirname(@__DIR__), "T1")))")
# Plots.scalefontsizes();
# ts_proportional_check()
# print_last_val()