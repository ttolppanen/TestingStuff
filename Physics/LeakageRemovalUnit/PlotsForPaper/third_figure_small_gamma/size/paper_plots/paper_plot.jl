using JLD2
using Plots

function t1_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T1", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "t", 
    ylabel = "First Site Population", 
    ylims = (0.2, 1)
    )
    for line in values(plotinfo.lines)
        size = parse(Float64, line.tag[4:end])
        plot!(pl, line.x, line.y, label = "Number of Transmons = " * string(size))
    end
    savefig(pl, joinpath(path, "T1.png"))
    savefig(pl, joinpath(path, "T1.svg"))
    return pl
end

function t2_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T2", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "t", 
    ylabel = "<σₓ> envelope", 
    ylims = (0.4, 1)
    )
    for line in values(plotinfo.lines)
        size = parse(Float64, line.tag[4:end])
        plot!(pl, line.x, line.y, label = "Number of Transmons = " * string(size))
    end
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
    return exp.(slope * collect(x) .+ c), -1 / slope 
end

function make_fit_line(dir, cut_size)
    path = joinpath(dirname(@__DIR__), dir)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")

    lines = values(plotinfo.lines)
    t_times = []
    t_time_x = []
    for line in lines
        fit_line_start = floor(Int, length(line.x) * cut_size)
        t_line, t_time = T1(line.x[fit_line_start:end], line.y[fit_line_start:end])

        t_x = parse(Float64, line.tag[5:end])
        push!(t_times, t_time)
        push!(t_time_x, t_x)
    end
    return t_time_x, t_times
end

function combine_fit_plot()
    t_time_xs, Ts = make_fit_line("Ts", 0.5)
    t_time_x1, T1 = make_fit_line("T1", 0.05)
    t_time_x2, T2 = make_fit_line("T2", 0.05)

    pl = plot(t_time_xs, Ts,
    xlabel = "L",
    ylabel = "T",
    label = "Ts",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9
    )
    plot!(pl, t_time_x1, T1, label = "T1")
    plot!(pl, t_time_x2, T2, label = "T2")
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))
    return pl
end

t1_plot()
t2_plot()
combine_fit_plot()