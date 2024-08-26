using PlotAndSave
using Plots
using JLD2

function T1(x, y) 
    # linreg for slope
    # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    res = (A \ log.(y))
    slope = res[1]
    c = res[2]
    return exp.(slope * collect(x) .+ c), -1 / slope 
end

function make_fit_plot(dir)
    path = joinpath(@__DIR__, dir)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")

    pl = plot(xlabel = plotinfo.xlabel,
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    title = plotinfo.title
    )

    lines = values(plotinfo.lines)
    t_times = []
    t_time_x = []
    for line in lines
        fit_line_start = floor(Int, length(line.x) * 0.6)
        t_line, t_time = T1(line.x[fit_line_start:end], line.y[fit_line_start:end])
        plot!(pl, line.x, line.y, label = line.tag)
        plot!(pl, line.x[fit_line_start:end], t_line, label = line.tag * " fit")

        t_x = parse(Float64, line.tag[5:end])
        push!(t_times, t_time)
        push!(t_time_x, t_x)
    end
    bar_x_label = collect(lines)[1].tag[1]
    bars = plot(t_time_x, t_times,
    xlabel = bar_x_label,
    ylabel = "Ts",
    dpi = 300,
    titlefontsize = 9,
    title = plotinfo.title,
    )
    savefig(pl, joinpath(path, "fit.png"))
    savefig(bars, joinpath(path, "T_times.png"))
    return pl
end

function make_fit_line(dir, cut_size)
    path = joinpath(@__DIR__, dir)
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
    t_time_x, Ts = make_fit_line("Ts", 0.5)
    t_time_x, T1 = make_fit_line("T1", 0.05)
    t_time_x, T2 = make_fit_line("T2", 0.05)

    pl = plot(t_time_x, Ts,
    xlabel = "d",
    ylabel = "T",
    label = "Ts",
    dpi = 300,
    titlefontsize = 9
    )
    plot!(pl, t_time_x, T1, label = "T1")
    plot!(pl, t_time_x, T2, label = "T2")
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
end

# include("Ts/code.jl")
# include("T1/code.jl")
# include("T2/code.jl")
# make_fit_plot("Ts/")
# make_fit_plot("T1/")
# make_fit_plot("T2/")
combine_fit_plot()
# exit()
