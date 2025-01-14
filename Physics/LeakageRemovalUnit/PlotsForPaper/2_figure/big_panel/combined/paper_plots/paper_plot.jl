using JLD2
using Plots
using PlotAndSave
include("/home/htolppan/Documents/Code/TestingStuff/Physics/LeakageRemovalUnit/PlotsForPaper/figure_colors.jl")

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([0.04, 50, 100, 10^-3, 1], [0.04, "50", "100", "10-3", "10^0"]),
    dpi = 300,
    xlabel = "Rate",
    ylabel = "Leakage Population in the Array", 
    ylims = (0, 1)
    )
    for line in values(plotinfo.lines)
        if line.tag == "gamma"
            tag = "Periodic Measurement"
            plot!(pl, line.x, line.y / 2, label = tag, color = color_periodic)
        elseif line.tag == "diss"
            tag = "Dissipation"
            plot!(pl, line.x, line.y / 2, label = tag, color = color_dissipation_removal)
        end
    end
    vline!(pl, [0.03], color = color_small_msr, ls = :dash, label = "")
    vline!(pl, [30], color = color_high_msr, ls = :dash, label = "")
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()