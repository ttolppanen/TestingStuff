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
        t_line, t_time = T1(line.x, line.y)
        plot!(pl, line.x, line.y, label = line.tag)
        plot!(pl, line.x, t_line, label = line.tag * " fit")

        t_x = parse(Float64, line.tag[5:end])
        push!(t_times, t_time)
        push!(t_time_x, t_x)
    end
    bar_x_label = collect(lines)[1].tag[1]
    bars = plot(t_time_x, t_times,
    seriestype = :bar,
    xlabel = bar_x_label,
    ylabel = "T",
    dpi = 300,
    titlefontsize = 9,
    title = plotinfo.title,
    )
    savefig(pl, joinpath(path, "fit.png"))
    savefig(bars, joinpath(path, "T_times.png"))
    return pl
end
include("gamma/small_gamma_U25/stack/code.jl")
include("gamma/small_gamma_U25/T1/code.jl")
include("gamma/small_gamma_U25/T2/code.jl")
make_fit_plot("gamma/small_gamma_U25/stack")
make_fit_plot("gamma/small_gamma_U25/T1")
make_fit_plot("gamma/small_gamma_U25/T2")