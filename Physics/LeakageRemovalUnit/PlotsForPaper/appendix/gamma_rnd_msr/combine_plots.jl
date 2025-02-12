using PlotAndSave
using Plots

function f()
    fig2_path = @__DIR__
    gamma_dt_001_path = joinpath(@__DIR__, "old_data", "gamma_dt_001")
    L2_data = joinpath(@__DIR__, "old_data", "L2_shorter_time")
    
    pl_fig2_line = loaddata(fig2_path).lines["L = 3"]
    pl_gamma_dt_001_line = loaddata(gamma_dt_001_path).lines["L = 3"]
    pl_L2_line = loaddata(L2_data).lines["L = 2"]

    pl = plot(xaxis = :log10)
    plot!(pl, pl_fig2_line.x, pl_fig2_line.y, label="Supplementary Fig. 2")
    # plot!(pl, pl_gamma_dt_001_line.x, pl_gamma_dt_001_line.y, label="L = 3 too small dt")
    plot!(pl, pl_L2_line.x, pl_L2_line.y, label="L = 2")

    vline!(pl, [50, 60], label = "50 and 60")
    vline!(pl, [0.08], label = "0.08")
    
    pl
end

f()