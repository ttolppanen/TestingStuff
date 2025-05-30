using Plots
using HDF5
using LaTeXStrings

function save_data_to_hdf5(plt, annonation)
    annonation = annonation[2]
    path = joinpath(@__DIR__, "hdf5_data", "fig3_data.h5")
    if annonation == 'a'
        file = h5open(path, "w")
    else
        file = h5open(path, "r+")
    end
    g = create_group(file, string(annonation))
    for s in plt.series_list
        label = String(s[:label])
        if label == "\$T_1\$"
            label = "t1"
        elseif label == "\$T_2\$"
            label = "t2"
        elseif label == "\$T_\\star\$" || label == "\$T_s\$"
            label = "ts"
        else
            continue
        end
        x = s[:x]
        y = s[:y]
        sg = create_group(g, label)
        write(sg, "y", y)
        write(sg, "x", x)
    end
    close(file)
end

function f(plots...)
    layout = @layout [
        a b c;
        d e f
    ]
    font_size = 10
    plots = deepcopy(collect(plots))
    annonations = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for i in 1:6
        save_data_to_hdf5(plots[i], annonations[i])
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