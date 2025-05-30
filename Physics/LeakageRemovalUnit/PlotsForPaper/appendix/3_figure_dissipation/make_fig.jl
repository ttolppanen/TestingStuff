using Plots
using HDF5
using LaTeXStrings

function save_data_to_hdf5(plt, annonation)
    annonation = annonation[2]
    path = joinpath(@__DIR__, "hdf5_data", "supplementary_fig3_data.h5")
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

function add_attributes_to_hdf5_data()
    path = joinpath(@__DIR__, "hdf5_data", "supplementary_fig3_data.h5")
    file = h5open(path, "r+")
    attrs(file)["DOI"] = "10.5281/zenodo.15516345"
    attrs(file)["description"] = "Supplementary figure 3 data for the low dissipation rate optimal spot. Groups contain the data of each subplot."
    attrs(file)["y_description"] = "Decay times T1, T2 and Ts in the units of microseconds."

    g_a = file["a"]
    g_b = file["b"]
    g_c = file["c"] 
    g_d = file["d"]
    g_e = file["e"]
    g_f = file["f"]

    attrs(g_a)["description"] = "T1, T2 and Ts as a function of disorder."
    attrs(g_a)["x_description"] = "Disorder over two pi in the units of MHz."

    attrs(g_b)["description"] = "T1, T2 and Ts as a function of coupling."
    attrs(g_b)["x_description"] = "Coupling over two pi in the units of MHz."

    attrs(g_c)["description"] = "T1, T2 and Ts as a function of array length."
    attrs(g_c)["x_description"] = "Number of transmons in the system, including the coding qubit."

    attrs(g_d)["description"] = "T1, T2 and Ts as a function of disorder in experimental conditions."
    attrs(g_d)["x_description"] = "Disorder over two pi in the units of MHz."

    attrs(g_e)["description"] = "T1, T2 and Ts as a function of coupling in experimental conditions."
    attrs(g_e)["x_description"] = "Coupling over two pi in the units of MHz."

    attrs(g_f)["description"] = "T1, T2 and Ts as a function of array length in experimental conditions."
    attrs(g_f)["x_description"] = "Number of transmons in the system, including the coding qubit."
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
    savefig(pl, joinpath(@__DIR__, "decay_times_dissipation_low.svg"))
    return pl
end

# include(joinpath(@__DIR__, "disorder/paper_plots/paper_plot.jl"))
# pl_diss = combine_fit_plot()

# include(joinpath(@__DIR__, "disorder_realistic/paper_plots/paper_plot.jl"))
# pl_diss_real = combine_fit_plot()

# include(joinpath(@__DIR__, "J/paper_plots/paper_plot.jl"))
# pl_J = combine_fit_plot()

# include(joinpath(@__DIR__, "J_realistic/paper_plots/paper_plot.jl"))
# pl_J_real = combine_fit_plot()

# include(joinpath(@__DIR__, "size/paper_plots/paper_plot.jl"))
# pl_L = combine_fit_plot()

# include(joinpath(@__DIR__, "size_realistic/paper_plots/paper_plot.jl"))
# pl_L_real = combine_fit_plot()

# f(pl_diss, pl_J, pl_L, pl_diss_real, pl_J_real, pl_L_real)
add_attributes_to_hdf5_data()