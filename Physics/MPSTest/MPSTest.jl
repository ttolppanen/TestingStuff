using MIPTDataProcessing
using Gadfly

function f()
    filename = "mipt_same_results"
    groups = ["mipt_15", "mipt_16"]
    for g in groups
        println(get_attributes_string(filename, g; addnewline = false))
    end
    layers = [getlayer(g) for g in groups]
    title = get_attributes_string(filename, groups[1])
    p = plot(layers..., Guide.title(title), Guide.xlabel("p"))
    p |> SVG("test.svg", 6inch, 4inch)
end

function getlayer(groupname)
    filename = "mipt_same_results"
    prob = get_probabilities(filename, groupname)
    dataend = get_data_end_mean(filename, groupname, prob, :entanglement_half)
    L = get_attribute(filename, groupname, "L")
    newLayer = layer(x=prob, y=dataend, Geom.line, color = ["L = $L"])
end

f()