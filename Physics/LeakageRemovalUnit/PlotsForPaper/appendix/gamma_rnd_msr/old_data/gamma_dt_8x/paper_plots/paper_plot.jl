using JLD2
using Plots
using PlotAndSave

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.03, 30], [0.03, "30"]),
    dpi = 300,
    xlabel = "Measurement Rate",
    ylabel = "Leakage Population in the Array", 
    # ylims = (10^-5, 1)
    )
    # vline!(pl, [0.009, 0.04, 3, 10])
    line = plotinfo.lines["L = 3"]
    plot!(pl, line.x, line.y, label = "")
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()