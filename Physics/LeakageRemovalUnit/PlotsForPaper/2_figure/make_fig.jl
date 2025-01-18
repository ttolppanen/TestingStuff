using Plots

function f(plot1 , plot2, plotbottom)
    layout = @layout [
        a b;
        c{0.66h}
    ]
    plot1 = plot(plot1, topmargin = 6.1Plots.mm)
    plot2 = plot(plot2, topmargin = 6.1Plots.mm)
    annotate!(plot1, ((-0.1, 1.22), text("(a)", 10, "Computer Modern")))
    annotate!(plot2, ((-0.1, 1.22), text("(b)", 10, "Computer Modern")))
    annotate!(plotbottom, ((-0.05, 1.15), text("(c)", 10, "Computer Modern")))
    pl = plot(plot1, plot2, plotbottom,
    layout = layout,
    size = (500, 500)
    )
    savefig(pl, joinpath(@__DIR__, "fig2.png"))
    savefig(pl, joinpath(@__DIR__, "fig2.svg"))
    return pl
end

include(joinpath(@__DIR__, "first_panel/time/paper_plots/paper_plot.jl"))
plot1 = make_plot()

include(joinpath(@__DIR__, "second_panel/time/paper_plots/paper_plot.jl"))
plot2 = make_plot()

include(joinpath(@__DIR__, "big_panel/combined/paper_plots/paper_plot.jl"))
plotbottom = make_plot()

f(plot1, plot2, plotbottom)
