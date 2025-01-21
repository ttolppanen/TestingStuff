using JLD2
using Plots
using PlotAndSave
using LRUFigureColors
using LaTeXStrings

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    xticks = ([10^-5, 10^-4, 10^-3, 10^0, 10^1, 10^2, 0.03, 30], ["10^{-5}", "10^{-4}", "10^{-3}", "10^{0}", "10^{1}",  "10^{2}", "0.03", "30"]),
    dpi = 300,
    xlabel = L"$\Gamma / J$",
    # ylabel = L"$\sum_{\ell} \left\langle \hat{n}_\ell(\hat{n}_\ell-\hat{I})/2 \right\rangle$",
    ylabel = L"$\textrm{Array \ Leakage} \ \left( P_2^{(\ell)} \right)$",
    legend = :bottomleft,
    ylims = (0, 1)
    )
    for line in values(plotinfo.lines)
        if line.tag == "gamma"
            tag = "Periodic Measurement"
            plot!(pl, line.x, line.y / 2, label = tag, color = color_periodic)
        elseif line.tag == "diss"
            tag = "Dissipation"
            plot!(pl, line.x, line.y / 2, label = tag, color = color_dissipation_removal, ls = :dash)
        end
    end
    vline!(pl, [0.03], color = color_small_msr, ls = :dash, label = "")
    vline!(pl, [30], color = color_high_msr, ls = :dash, label = "")
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()