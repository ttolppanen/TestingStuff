using JLD2
using Plots

function make_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "t", 
    ylabel = "First Site Leakage Population", 
    # ylims = (10^-5, 1)
    )
    for line in values(plotinfo.lines)
        if line.tag != "gamma = 0.1"
            gamma = parse(Float64, line.tag[8:end])
            if gamma != 0.03
                gamma = Int(floor(gamma))
            end
            plot!(pl, line.x, line.y, label = "Measurement Rate = " * string(gamma))
        end
    end
    savefig(pl, joinpath(path, "paper_plot.png"))
    savefig(pl, joinpath(path, "paper_plot.svg"))
    return pl
end

make_plot()