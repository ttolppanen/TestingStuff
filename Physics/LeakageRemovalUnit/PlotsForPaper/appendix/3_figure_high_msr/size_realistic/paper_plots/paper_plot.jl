using JLD2
using Plots
using OrderedCollections
using PlotAndSave
using QuantumTimes
using LaTeXStrings
using LRUFigureColors


function TS_fit()
    path = joinpath(dirname(@__DIR__), "Ts")
    T_fit(path, 0.001)
end
function T1_fit()
    path = joinpath(dirname(@__DIR__), "T1")
    T_fit(path, 0.05)
end
function T2_fit()
    path = joinpath(dirname(@__DIR__), "T2")
    T_fit(path, 0.05, 0.9)
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
        size = parse(Int64, line.tag[4:end])
        plot!(pl, line.x, line.y, label = "L = " * string(size))
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
        size = parse(Int64, line.tag[4:end])
        plot!(pl, line.x, line.y, label = "L = " * string(size))
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
    ylabel = L"$\langle \sigma_x \rangle \textrm{envelope}$"
    )
    fit_lines = T2_fit()
    for line in values(plotinfo.lines)
        size = parse(Int64, line.tag[4:end])
        plot!(pl, line.x, line.y, label = "L = " * string(size))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    savefig(pl, joinpath(path, "T2.png"))
    savefig(pl, joinpath(path, "T2.svg"))
    return pl
end

function T1(x, y) 
    # linreg for slope
    # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    res = (A \ log.(y))
    slope = res[1]
    c = res[2]
    return x -> exp.(slope * x .+ c), -1 / slope 
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
    t_time_xs = [parse(Int64, key[4:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    lines = T1_fit()
    t_time_x1 = [parse(Int64, key[4:end]) for key in keys(lines)]
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    t_time_x2 = [parse(Int64, key[4:end]) for key in keys(lines)]
    T2 = [val[2] for val in values(lines)]

    path = joinpath(dirname(@__DIR__), "T1")
    J = get_param(path, :J) * 2 * pi
    pl = plot(t_time_xs, Ts / J, seriestype = :scatter,
    xlabel = L"$L$",
    ylabel = L"$\textrm{Survival\ Times}\ [\mu \textrm{s}]$",
    label = L"$T_s$",
    color = color_ts,
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9,
    ylims = [0, 17],
    yticks = 0:4:16,
    yminorticks = 4,
    xminorticks = 1,
    minorgrid = true
    )
    plot!(pl, t_time_x1, T1 / J, seriestype = :scatter, label = L"$T_1$", color = color_t1)
    plot!(pl, t_time_x2, T2 / J, seriestype = :scatter, label = L"$T_2$", color = color_t2)

    
    t1 = get_param(path, :diss)
    hline!(pl, [t1], ls = :dash, label =  L"$T_1^{\textrm{qb}}$", color = color_t1)
    pure_t2 = get_param(path, :deco)
    t2_real = t2(t1, pure_t2)
    hline!(pl, [t2_real], ls = :dash, label =  L"$T_2^{\textrm{qb}}$", color = color_t2)

    U = get_param(path, :U) * 2 * pi
    short_scale(L) = t_small_scale(J, U, L)
    t_y = short_scale.(t_time_xs)
    plot!(pl, t_time_xs, t_y, label = "Short", ls = :dash, color = color_ts_short)

    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))
    return pl
end

Plots.scalefontsizes()
ts_plot()
t1_plot()
t2_plot()
combine_fit_plot();
println("Trajectories: $(get_traj(joinpath(dirname(@__DIR__), "T1")))")
Plots.scalefontsizes()