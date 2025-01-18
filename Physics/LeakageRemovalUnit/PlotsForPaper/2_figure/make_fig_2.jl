function f(plot1 , plot2, plotbottom)
    layout = @layout [
        a b;
        c{0.66h}
    ]
    plot1 = plot(plot1, topmargin = 7Plots.mm)
    plot2 = plot(plot2, topmargin = 7Plots.mm)
    annotate!(plot1, ((-0.08, 1.22), "a)"))
    annotate!(plot2, ((-0.08, 1.22), "b)"))
    annotate!(plotbottom, ((-0.04, 1.15), "c)"))
    pl = plot(plot1, plot2, plotbottom,
    layout = layout,
    size = (500, 500)
    )
    path = dirname(@__FILE__)
    savefig(pl, joinpath(path, "fig2.png"))
    savefig(pl, joinpath(path, "fig2.svg"))
    return pl
end

include(joinpath(@__DIR__, "first_panel/time/paper_plots/paper_plot.jl"))
plot1 = make_plot()

include(joinpath(@__DIR__, "second_panel/time/paper_plots/paper_plot.jl"))
plot2 = make_plot()

include(joinpath(@__DIR__, "big_panel/combined/paper_plots/paper_plot.jl"))
plotbottom = make_plot()

f(plot1, plot2, plotbottom)
