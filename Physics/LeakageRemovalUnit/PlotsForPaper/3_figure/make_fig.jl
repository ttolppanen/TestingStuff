using Plots

function f(plots...)
    layout = @layout [
        a b c;
        d e f
    ]
    font_size = 10
    plots = deepcopy(collect(plots))
    # pl_diss = plot(plot_diss, topmargin = 6.1Plots.mm)
    # plot2 = plot(plot2, topmargin = 6.1Plots.mm)
    annonations = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for i in 1:6
        plot!(plots[i], topmargin = 7Plots.mm, ytickfont = font_size - 1, xtickfont = font_size - 1)
        annotate!(plots[i], ((-0.01, 1.15), text(annonations[i], font_size, "Computer Modern")))

        if !(i == 1 || i == 4)
            plot!(plots[i], legend = :none)
        end
    end
    
    pl = plot(plots...,
    layout = layout,
    guidefont = font_size,
    legendfont = 7,
    margin = 3Plots.mm,
    size = (1000, 500)
    )
    savefig(pl, joinpath(@__DIR__, "fig3.png"))
    savefig(pl, joinpath(@__DIR__, "fig3.svg"))
    return pl
end

include(joinpath(@__DIR__, "disorder/paper_plots/paper_plot.jl"))
pl_diss = combine_fit_plot()

include(joinpath(@__DIR__, "disorder_realistic/paper_plots/paper_plot.jl"))
pl_diss_real = combine_fit_plot()

include(joinpath(@__DIR__, "J/paper_plots/paper_plot.jl"))
pl_J = combine_fit_plot()

include(joinpath(@__DIR__, "J_realistic/paper_plots/paper_plot.jl"))
pl_J_real = combine_fit_plot()

include(joinpath(@__DIR__, "size/paper_plots/paper_plot.jl"))
pl_L = combine_fit_plot()

include(joinpath(@__DIR__, "size_realistic/paper_plots/paper_plot.jl"))
pl_L_real = combine_fit_plot()

f(pl_diss, pl_J, pl_L, pl_diss_real, pl_J_real, pl_L_real)