using PlotAndSave
using Gadfly
using JLD2

function plot_to_log(name)
    path = joinpath(@__DIR__, name)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    layers = []
    for line in values(plotinfo.lines)
        log_layer = layer(x=line.x, y=line.y, Geom.line, color = [line.tag])
        push!(layers, log_layer)
    end
    path = joinpath(path, "log plot.svg")
    p = plot(layers..., Scale.y_log10, Guide.title(plotinfo.title), Guide.xlabel(plotinfo.xlabel), Guide.ylabel(plotinfo.ylabel))
    p |> SVG(path, 6inch, 4inch)
end

plot_to_log("mean_steps_scramble0.1_even larger_systems")