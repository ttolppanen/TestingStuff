using Plots
using PlotAndSave

function make_loss_comparison(path)
    println(path)
    plotinfo = loaddata(path)

    T1_best = 1 / plotinfo.parameters[:diss]
    T2_best = 2 / plotinfo.parameters[:deco]

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
end

make_loss_comparison(joinpath(dirname(@__DIR__), "paper_plots", "data.jld2"))