using PlotAndSave
using JLD2
using Plots

function f()
    path = @__DIR__
    plotinfo_normal = load(joinpath(path, "data.jld2"), "plotinfo")
    plotinfo_large = load(joinpath(path, "large_diss_zeno", "data.jld2"), "plotinfo")
    x = vcat(plotinfo_normal.lines["L = 3"].x, plotinfo_large.lines["L = 3"].x[2:end])
    y = vcat(plotinfo_normal.lines["L = 3"].y, plotinfo_large.lines["L = 3"].y[2:end])
    plot(x, y, xaxis = :log10, dpi = 300,
    xticks = ([0.01, 0.03, 0.1, 5, 10, 20, 50, 100, 200], [0.01, 0.03, 0.1, 5, 10, 20, 50, 100, 200]),)
end

f()