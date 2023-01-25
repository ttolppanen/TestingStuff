include("./PAPER_ReadData.jl")
using MIPTDataProcessing

function f()
    oldlayer = MakePlotFromOldData()
    filename = "mipt_same_results"
    groupname = "mipt_5"
    prob = get_probabilities(filename, groupname)
    dataend = get_data_end_mean(filename, groupname, prob, :entanglement_half)
    newLayer = layer(x=prob, y=dataend, Geom.line, color = ["New Data"])
    title = get_attributes_string(filename, groupname)
    cut1 = 61
    cut2 = cut1 + 21
    title = title[1:cut1] * "\n" * title[cut1:cut2] * "\n" * title[cut2:end]
    p = plot(oldlayer, newLayer, Guide.title(title))
    p |> SVG("test.svg", 6inch, 4inch)
end

f()