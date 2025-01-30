using PlotAndSave
using Plots

function f()
    path1 = joinpath(@__DIR__, "gamma_rnd_msr_smaller_dt/")
    path2 = joinpath(@__DIR__, "gamma_too_large_dt/")
    new_path = joinpath(@__DIR__, "dt_check")

    plot1 = loaddata(path1)
    plot2 = loaddata(path2)
    plot1.lines["L = 3"]

    pl = plot()
    line = plot1.lines["L = 3"]
    plot!(pl, line.x, line.y, xaxis = :log10, label = "dt = 0.005")
    line = plot2.lines["L = 3"]
    plot!(pl, line.x, line.y, xaxis = :log10, label = "dt = 0.05")
    vline!(pl, [0.055, 50])
end

f()