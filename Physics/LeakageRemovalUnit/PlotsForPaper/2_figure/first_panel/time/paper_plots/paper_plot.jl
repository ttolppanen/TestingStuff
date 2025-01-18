using JLD2
using Plots
using LRUFigureColors
using PlotAndSave
using LaTeXStrings

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = L"$\textrm{Time} \ \left( \frac{t}{J} \right)$", 
    ylabel = L"$\left\langle \hat{n}_1(\hat{n}_1-\hat{I})/2 \right\rangle$"
    # ylims = (10^-5, 1)
    )
    l_low = plotinfo.lines["gamma = 0.03"]
    l_high = plotinfo.lines["gamma = 50.0"]
    plot!(pl, l_low.x, l_low.y / 2, label = L"$ \Gamma_{\textrm{fb}} = 0.03J$", color = color_small_msr)
    plot!(pl, l_high.x, l_high.y / 2, label = L"$\Gamma_{\textrm{fb}}  = 50J$", color = color_high_msr)
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()