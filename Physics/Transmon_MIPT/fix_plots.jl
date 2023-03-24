using PlotAndSave
using JLD2
using Gadfly

function f(filename)
    path = joinpath(@__DIR__, filename, "data.jld2")
    plotinfo = load(path, "plotinfo")
    lines = 
    [
        plotinfo.lines["maxdim = 30"], 
        plotinfo.lines["maxdim = 50"], 
        plotinfo.lines["maxdim = 70"], 
        plotinfo.lines["maxdim = 90"], 
        plotinfo.lines["maxdim = 110"]
    ]
    lines = reverse(lines)
    path = joinpath(@__DIR__, "fixed_" * filename * ".svg")
    layers = []
    for lineinfo in lines
        new_layer = layer(x=lineinfo.x, y=lineinfo.y, Geom.line, color = [lineinfo.tag])
        push!(layers, new_layer)
    end
    p = plot(layers..., Guide.title(""), Guide.xlabel(plotinfo.xlabel), Guide.ylabel(plotinfo.ylabel))
    p |> SVG(path, 5.6inch, 3.3inch)
end

f("maxdim_fluc")