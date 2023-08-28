using PlotAndSave
using JLD2
using Plots
using Statistics

function average_fluctuations()
    path = joinpath(@__DIR__, "fluc_test")
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = "p", 
    ylabel = plotinfo.ylabel, 
    dpi = 300, 
    legend = :outertopright,
    titlefontsize = 9
    )
    x = []
    y = []
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        push!(x, p)
        x_ticks = length(line.y)
        push!(y, mean(line.y[Int(floor(x_ticks / 2)) : end]))
    end
    plot!(pl, x, y)
    path = joinpath(path, "average fluctuations.png")
    savefig(pl, path)
    return pl
end

average_fluctuations()