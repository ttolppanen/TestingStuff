using JLD2
using Plots
using PlotAndSave
using LRUFigureColors
using LaTeXStrings

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    prd_plotinfo = loaddata(path, "periodic_data.jld2")
    pl = plot(xaxis = :log10,
    xticks = ([10^-5, 10^-4, 10^-3, 10^-2, 10^0, 10^1, 0.055, 60], ["10^{-5}", "10^{-4}", "10^{-3}", "10^{-2}", "10^{0}", "10^{1}", "0.055", "60"]),
    dpi = 300,
    xlabel = L"$\Gamma / J$",
    # ylabel = L"$\sum_{\ell} \left\langle \hat{n}_\ell(\hat{n}_\ell-\hat{I})/2 \right\rangle$",
    ylabel = L"$\textrm{Array \ Leakage} \ P_\star^{(L)}$",
    legend = :bottomleft,
    legendfontsize = 7,
    ylims = (0, 1)
    )
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y / 2, label = "Random feedback measurement", color = color_ts_short)
    end
    line = prd_plotinfo.lines["L = 3"]
    plot!(pl, line.x, line.y / 2, label = "Periodic feedback measurement", ls = :dash, color = color_periodic)
    vline!(pl, [0.055], color = color_small_msr, ls = :dash, label = "")
    vline!(pl, [60], color = color_high_msr, ls = :dash, label = "")

    annotate!(pl, ((0.55, 0.9), text(L"\Gamma_{\textrm{rnd}}^{\textrm{low}}", 10, "Computer Modern")))
    annotate!(pl, ((0.943, 0.9), text(L"\Gamma_{\textrm{rnd}}^{\textrm{high}}", 10, "Computer Modern")))

    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "rnd_msr_optimal_spots.svg"))
    return pl
end

make_plot()