using PlotAndSave
using Plots
using LRUFigureColors

function pl_path(name)
    return joinpath(dirname(dirname(@__DIR__)), "3_figure", name, "paper_plots", "paper_plot.jl")
end
function every_n_plot(pl, n)
    total_series = length(pl.series_list) - 1
    step = max(1, div(total_series, n))
    new_series_list = []
    for i in 1:step:total_series
        push!(new_series_list, pl.series_list[i])
        push!(new_series_list, pl.series_list[i + 1])
    end
    new_plot = plot()
    for series in new_series_list
        plot!(new_plot, series)
    end
    return new_plot
end

function remove_every_n!(plots, n)
    for i in eachindex(plots)
        plots[i] = every_n_plot(plots[i], n)
    end
end

function format_ts(ts, xlims)
    plot!(ts, yformatter = x -> x/2, xlims = xlims) 
end
function format_t1(t1, xlims)
    plot!(t1, xlims = xlims)
end
function format_t2(t2, xlims)
    cut_last_10_percent!(t2)
    plot!(t2, ylims = [0, 1], xlims = xlims)
end
function cut_last_10_percent!(plot)
    for series in plot.series_list
        n = length(series[:x])
        cut_index = floor(Int, 0.9 * n)
        series[:x] = series[:x][1:cut_index]
        series[:y] = series[:y][1:cut_index]
    end
end

function make_plots(ts, t1, t2, xlims)
    plots = [ts, t1, t2]
    remove_every_n!(plots, 4)
    format_ts(plots[1], xlims)
    format_t1(plots[2], xlims)
    format_t2(plots[3], xlims)
    for i in eachindex(plots)
        plot!(plots[i], legend = false)
    end
    pl = plot(plots..., layout = @layout([a; b; c]))
    return pl
end

function combine_plots(plots)
    combined_plot = plot(plots..., 
    layout = @layout([a b c; d e f]),
    size = (2000, 2000)
    )
    return combined_plot
end

plots = []

include(pl_path("disorder"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 300]))

include(pl_path("J"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 10000]))

include(pl_path("size"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 300]))

include(pl_path("disorder_realistic"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 300]))

include(pl_path("J_realistic"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 10000]))

include(pl_path("size_realistic"))
push!(plots, make_plots(ts_plot(), t1_plot(), t2_plot(), [0, 300]))

combine_plots(plots)
