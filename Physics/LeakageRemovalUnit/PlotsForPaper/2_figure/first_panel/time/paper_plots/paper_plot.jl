using JLD2
using Plots
include("/home/htolppan/Documents/Code/TestingStuff/Physics/LeakageRemovalUnit/PlotsForPaper/figure_colors.jl")

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "First Site Leakage Population"
    # ylims = (10^-5, 1)
    )
    Plots.scalefontsizes(1.5)
    l_low = plotinfo.lines["gamma = 0.03"]
    l_high = plotinfo.lines["gamma = 50.0"]
    plot!(pl, l_low.x, l_low.y / 2, label = "Measurement Rate = 0.03", color = color_small_msr)
    plot!(pl, l_high.x, l_high.y / 2, label = "Measurement Rate = 50", color = color_high_msr)
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    Plots.scalefontsizes()
    return pl
end

make_plot()