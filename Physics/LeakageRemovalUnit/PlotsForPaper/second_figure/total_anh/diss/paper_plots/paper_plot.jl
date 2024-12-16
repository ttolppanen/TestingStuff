using JLD2
using Plots
using PlotAndSave

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.04, 50, 100], [0.04, "50", "100"]),
    dpi = 300,
    xlabel = "Dissipation Rate",
    ylabel = "Leakage Population in the Array", 
    ylims = (0, 2)
    )
    line = plotinfo.lines["L = 3"]
    plot!(pl, line.x, line.y, label = "")
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()