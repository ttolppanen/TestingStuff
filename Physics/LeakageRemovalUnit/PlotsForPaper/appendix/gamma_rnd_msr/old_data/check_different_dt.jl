using PlotAndSave
using Plots

function f()
    path1 = joinpath(@__DIR__, "gamma_dt_4x/")
    path2 = joinpath(@__DIR__, "gamma_dt_8x/")
    new_path = joinpath(@__DIR__, "dt_check")

    plot1 = loaddata(path1)
    plot2 = loaddata(path2)
    plot1.lines["L = 3"]

    pl = plot()
    line = plot1.lines["L = 3"]
    plot!(pl, line.x, line.y, xaxis = :log10, label = "...")
    line = plot2.lines["L = 3"]
    plot!(pl, line.x, line.y, xaxis = :log10, label = "comparison")
    vline!(pl, [0.055, 50])
end

f()