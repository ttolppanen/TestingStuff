using PlotAndSave
using Plots

function f()
    L2 = loaddata(joinpath(@__DIR__, "L2_data.jld2"))
    L3 = loaddata(joinpath(@__DIR__, "L3_data.jld2"))
    
    pl = plot(
        L2.lines["T1 loss %"].x,
        L2.lines["T1 loss %"].y,
        label = "L=2, T1 loss %",
        xlims = (0.0, 1.0),
        xlabel = L2.xlabel
    )
    plot!(pl, L3.lines["T1 loss %"].x,
    L3.lines["T1 loss %"].y,
    label = "L=3, T1 loss %"
    )
    plot!(pl, L2.lines["T2 loss %"].x,
    L2.lines["T2 loss %"].y,
    label = "L=2, T2 loss %"
    )
    plot!(pl, L3.lines["T2 loss %"].x,
    L3.lines["T2 loss %"].y,
    label = "L=3, T2 loss %"
    )
end

f()