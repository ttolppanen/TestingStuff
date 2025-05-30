using Plots
using LaTeXStrings
using HDF5

function save_data_to_hdf5(plt, annonation)
    path = joinpath(@__DIR__, "hdf5_data", "fig2_data.h5")
    if annonation == 'a'
        file = h5open(path, "w")
    else
        file = h5open(path, "r+")
    end
    g = create_group(file, string(annonation))
    println(annonation)
    for s in plt.series_list
        label = String(s[:label])
        println(label)
        if label == "\$ \\Gamma_{\\textrm{fb}}^{\\textrm{low}} = 0.03J\$" || label == "\$\\Gamma_{\\textrm{fb}}^{\\textrm{low}} = 0.03J\$"
            label = "low measurement rate"
        elseif label == "\$\\Gamma_{\\textrm{fb}}^{\\textrm{high}}  = 30J\$" || label == "\$\\Gamma_{\\textrm{fb}}^{\\textrm{high}} = 30J\$"
            label = "high measurement rate"
        elseif label == "Feedback measurement"
            label = "feedback measurement"
        elseif label == "Dissipation"
            label = "dissipation"
        else
            continue
        end
        x = collect(s[:x])
        y = s[:y]
        sg = create_group(g, label)
        write(sg, "y", y)
        write(sg, "x", x)
    end
    close(file)
end

function add_attributes_to_hdf5_data()
    path = joinpath(@__DIR__, "hdf5_data", "fig2_data.h5")
    file = h5open(path, "r+")
    attrs(file)["DOI"] = "10.5281/zenodo.15516345"
    attrs(file)["description"] = "Figure 2 data. Groups contain the data of each subplot."
    
    g_a = file["a"]
    g_b = file["b"]
    g_c = file["c"]

    attrs(g_a)["description"] = "Qubit leakage population as a function of time. Time is in the units of one over coupling 1/J."
    attrs(g_b)["description"] = "Array leakage population as a function of time. Time is in the units of one over coupling 1/J."
    attrs(g_c)["description"] = "Array leakage population as a function of measurement/dissipation. Rate is in the units of the coupling J."
    close(file)
end

function f(plot1 , plot2, plotbottom)
    layout = @layout [
        a b;
        c{0.66h}
    ]
    plot1 = plot(plot1, topmargin = 6.3Plots.mm)
    plot2 = plot(plot2, topmargin = 6.3Plots.mm)
    annotate!(plot1, ((-0.01, 1.17), text(L"$\textbf{(a)}$", 10, "Computer Modern")))
    annotate!(plot2, ((-0.01, 1.17), text(L"$\textbf{(b)}$", 10, "Computer Modern")))
    annotate!(plotbottom, ((-0.002, 1 + 0.1), text(L"$\textbf{(c)}$", 10, "Computer Modern")))

    annotate!(plotbottom, ((0.53, 0.9), text(L"\Gamma_{\textrm{fb}}^{\textrm{low}}", 10, "Computer Modern")))
    annotate!(plotbottom, ((0.93, 0.9), text(L"\Gamma_{\textrm{fb}}^{\textrm{high}}", 10, "Computer Modern")))

    plot!(plotbottom, topmargin = -1.5Plots.mm)
    pl = plot(plot1, plot2, plotbottom,
    layout = layout,
    size = (500, 500)
    )
    savefig(pl, joinpath(@__DIR__, "fig2.png"))
    savefig(pl, joinpath(@__DIR__, "fig2.svg"))
    for i in 1:3
        plt = [plot1, plot2, plotbottom][i]
        annotation = ['a', 'b', 'c'][i]
        save_data_to_hdf5(plt, annotation)
    end
    add_attributes_to_hdf5_data()
    return pl
end

include(joinpath(@__DIR__, "first_panel/time/paper_plots/paper_plot.jl"))
plot1 = make_plot()

include(joinpath(@__DIR__, "second_panel/time/paper_plots/paper_plot.jl"))
plot2 = make_plot()

include(joinpath(@__DIR__, "big_panel/combined/paper_plots/paper_plot.jl"))
plotbottom = make_plot()

f(plot1, plot2, plotbottom)
