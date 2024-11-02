using Plots
using PlotAndSave

function make_loss_comparison(path)
    println(path)
    plotinfo = loaddata(path)

    to_flt(x) = parse(Float64, match(r"\d+", x).match)
    U = (to_flt(plotinfo.parameters[:U]))
    T1_best = 1 / (to_flt(plotinfo.parameters[:diss]) * 10^(-3) / U)
    T2_best = 2 / (to_flt(plotinfo.parameters[:deco]) * 10^(-3) / U)

    Ts = plotinfo.lines["Ts"]
    T1 = plotinfo.lines["T1"]
    T2 = plotinfo.lines["T2"]

    T1_loss = 1 .- T1.y ./ T1_best
    T2_loss = 1 .- T2.y ./ T2_best
    Ts_perf = Ts.y ./ T2.y

    pl = plot(
        xlabel = plotinfo.xlabel,
        ylimits = (0, 1)
    )
    x = plotinfo.lines["Ts"].x
    plot!(pl, x, T1_loss, label = "T1 loss %")
    plot!(pl, x, T2_loss, label = "T2 loss %")
    plot!(pl, x, Ts_perf, label = "Ts / T2")

    savefig(pl, joinpath(@__DIR__, "plot.png"))
    savefig(pl, joinpath(@__DIR__, "plot.svg"))
end

make_loss_comparison(joinpath(dirname(@__DIR__), "paper_plots", "data.jld2"))