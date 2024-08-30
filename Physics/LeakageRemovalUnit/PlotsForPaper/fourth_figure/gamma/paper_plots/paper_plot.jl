using JLD2
using Plots

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.03, 45], [0.03, 45]),
    dpi = 300,
    xlabel = "Measurement Rate",
    ylabel = "Leakage Population in the Array", 
    # ylims = (10^-5, 1)
    )
    line = plotinfo.lines["L = 5"]
    plot!(pl, line.x, line.y, label = "")
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()