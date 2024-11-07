using JLD2
using Plots
using OrderedCollections

function TS_fit()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "Ts", "data.jld2"), "plotinfo")
    lines = OrderedDict()

    line = plotinfo.lines["dw = 5"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 10"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 15"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 20"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)
    
    return lines
end

function T1_fit()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T1", "data.jld2"), "plotinfo")
    lines = OrderedDict()

    line = plotinfo.lines["dw = 5"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 10"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 15"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 20"]
    f_s = floor(Int, length(line.x) * 0.01)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)
    
    return lines
end

function T2_fit()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T2", "data.jld2"), "plotinfo")
    lines = OrderedDict()

    line = plotinfo.lines["dw = 5"]
    f_s = floor(Int, length(line.x) * 0.25)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 10"]
    f_s = floor(Int, length(line.x) * 0.25)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 15"]
    f_s = floor(Int, length(line.x) * 0.25)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)

    line = plotinfo.lines["dw = 20"]
    f_s = floor(Int, length(line.x) * 0.25)
    f_e = floor(Int, length(line.x) * 1.0)
    x = line.x[f_s:f_e]
    y = line.y[f_s:f_e]
    lines[line.tag] = T1(x, y)
    
    return lines
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
        dw = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(dw))
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
        dw = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(dw))
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
    ylabel = "<σₓ> envelope"
    )
    fit_lines = T2_fit()
    for line in values(plotinfo.lines)
        dw = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(dw))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end

    # normal decoherence
    t = plotinfo.lines["dw = 5"].x
    deph_rate = plotinfo.parameters[:deco]
    plot!(pl, t, x -> exp(-x * deph_rate * 0.5), label = "decoherence")

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
    t_time_xs = [parse(Int64, key[5:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    lines = T1_fit()
    t_time_x1 = [parse(Int64, key[5:end]) for key in keys(lines)]
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    t_time_x2 = [parse(Int64, key[5:end]) for key in keys(lines)]
    T2 = [val[2] for val in values(lines)]

    pl = plot(t_time_xs, Ts,
    xlabel = "Disorder",
    ylabel = "Survival Times",
    label = "Tₛ",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9
    )
    plot!(pl, t_time_x1, T1, label = "T₁")
    plot!(pl, t_time_x2, T2, label = "T₂")
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))

    J = 5 * 2 * pi #5MHz
    Ts_best = Ts[end] / J # microseconds
    T1_best = T1[end] / J
    T2_best = T2[end] / J
    println("Ts = $(round(Ts_best, digits = 2))μs")
    println("T1 = $(round(T1_best, digits = 2))μs")
    println("T2 = $(round(T2_best, digits = 2))μs")
    return pl
end

Plots.scalefontsizes(1.5)
ts_plot()
t1_plot()
t2_plot()
combine_fit_plot();
Plots.scalefontsizes();