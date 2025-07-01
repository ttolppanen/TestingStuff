using PlotAndSave
using Plots
using LRUFigureColors
using LaTeXStrings

function pl_path(name)
    return joinpath(dirname(dirname(@__DIR__)), "3_figure", name, "paper_plots", "paper_plot.jl")
end

function last_two_plots(pl)
    new_series_list = pl.series_list[end-2:end-1]
    new_plot = plot()
    for series in new_series_list
        plot!(new_plot, series)
    end
    return new_plot
end

function format_ts(ts, xlims)
    plot!(ts, ylims = [0, 2.2], yformatter = x -> x/2, xlims = xlims) 
end

function format_t1(t1, xlims)
    plot!(t1, ylims = [0, 1], xlims = xlims)
end

function format_t2(t2, xlims)
    # <cut_last_10_percent!(t2[end])
    plot!(t2, ylims = [0, 1], xlims = xlims)
end

function cut_last_10_percent!(plot)
    series = plot.series_list[1]
    n = length(series[:x])
    start_cut = 1 + floor(Int, 0.002 * n)
    end_cut = floor(Int, 0.9 * n)
    series[:x] = series[:x][start_cut:end_cut]
    series[:y] = series[:y][start_cut:end_cut]
end

# Ts: The label is "$T_\star$".
# T1: The label is "$T_1$".
# T2: The label is "$T_2$".

function plots_to_lines(t1, t2, ts)
    new_plot = plot()
    t1.series_list[1][:label] = L"\langle n_1 \rangle"
    t1.series_list[1][:linecolor] = color_t1
    plot!(new_plot, t1.series_list[1])
    plot!(new_plot, t1.series_list[2])
    t2.series_list[1][:label] = L"\langle \sigma_x \rangle"
    t2.series_list[1][:linecolor] = color_t2
    plot!(new_plot, t2.series_list[1])
    plot!(new_plot, t2.series_list[2])
    ts.series_list[1][:label] = L"P_\star^{(L)}"
    ts.series_list[1][:linecolor] = color_ts
    ts.series_list[1][:y] *= 0.5
    ts.series_list[2][:y] *= 0.5
    plot!(new_plot, ts.series_list[1])
    plot!(new_plot, ts.series_list[2])
    return new_plot
end

function make_plots(t1, t2, ts, title)
    plots = [t1, t2, ts]
    for i in eachindex(plots)
        plots[i] = last_two_plots(plots[i])
    end
    for i in eachindex(plots)
        plot!(plots[i], legend = false)
    end
    cut_last_10_percent!(plots[2])
    pl = plots_to_lines(plots[1], plots[2], plots[3])
    plot!(pl, title = title, xlabel = L"tJ")
    
    if title == "Ideal"
        plot!(pl, legend = :right)
    end

    t = t_small_scale(1, 50) * 2 # U / J = 50, L - 1 = 2
    vline!(pl, [t], ls = :dash, color = color_ts_short_fig2, label = L"(L-1)T_{\textrm{prop}}")
    return pl
end

function combine_plots(plots)
    combined_plot = plot(plots..., 
    layout = @layout([a b]),
    size = (800, 400) .* 0.7
    )
    return combined_plot
end

plots = []

include(pl_path("disorder"))
push!(plots, make_plots(t1_plot(), t2_plot(), ts_plot(), "Ideal"))

include(pl_path("disorder_realistic"))
push!(plots, make_plots(t1_plot(), t2_plot(), ts_plot(), "Experimental Conditions"))

pl = combine_plots(plots)
savefig(joinpath(@__DIR__, "fit.png"))
savefig(joinpath(@__DIR__, "fig3_fits.svg"))
pl
