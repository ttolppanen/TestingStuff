using Plots

function f(plots...)
    layout = @layout [
        a b c;
        d e f
    ]
    font_size = 10
    plots = deepcopy(collect(plots))
    annonations = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for i in 1:6
        plot!(plots[i], topmargin = 7Plots.mm, ytickfont = font_size - 1, xtickfont = font_size - 1)
        annotate!(plots[i], ((-0.00, 1.15), text(latexstring("\\textbf{$(annonations[i])}"), font_size, "Computer Modern")))

        if (i >= 1 && i <= 3)
            plot!(plots[i], bottommargin = -3Plots.mm)
        end
        if !(i == 1 || i == 4)
            plot!(plots[i], legend = :none, ylabel = "")
            plot!(plots[i], leftmargin = -3Plots.mm)
        end
        if i == 4
            plot!(plots[i], 
                background_color_legend = RGBA(1, 1, 1, 0.66),
                foreground_color_legend = RGBA(0, 0, 0, 1)
            )
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
    savefig(pl, joinpath(@__DIR__, "decay_times_dissipation_high.svg"))
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