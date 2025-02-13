using PlotAndSave
using Plots

function f()
    fig2_path = @__DIR__
    gamma_dt_001_path = joinpath(@__DIR__, "old_data", "rnd_msr_no_dissorder")
    
    pl_fig2_line = loaddata(fig2_path).lines["L = 3"]
    pl_gamma_dt_001_line = loaddata(gamma_dt_001_path).lines["L = 3"]

    pl = plot(xaxis = :log10)
    plot!(pl, pl_fig2_line.x, pl_fig2_line.y, label="Supplementary Fig. 2")
    plot!(pl, pl_gamma_dt_001_line.x, pl_gamma_dt_001_line.y, label="L = 3 no dissorder")

    # vline!(pl, [50, 60], label = "50 and 60")
    # vline!(pl, [0.08], label = "0.08")
    
    pl
end

f()