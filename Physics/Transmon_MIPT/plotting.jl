using MIPTDataProcessing
using Gadfly

function f()
    filename = "mipt_same_results"
    # groups = ["mipt_33", "mipt_34", "mipt_35"]
    groups = ["mipt_36", "mipt_37", "mipt_38"]
    for g in groups
        println(get_attributes_string(filename, g; addnewline = false))
    end
    dtype = "entanglement_half"
    layers = [getlayer(g, dtype) for g in groups]
    title = get_attributes_string(filename, groups[1])
    p = plot(layers..., Guide.title(title), Guide.xlabel("p"))
    p |> SVG(dtype * "_40.svg")
end

function getlayer(groupname, dataset)
    filename = "mipt_same_results"
    dt = get_attribute(filename, groupname, "dt")
    prob = get_probabilities_dataset(filename, groupname, dataset)
    # prob = filter(p -> p > 0 && p < 0.05, prob)
    # prob = filter(p -> p / dt > 0.49 && p / dt < 2.01, prob)
    dataend = get_data_end_mean(filename, groupname, prob, dataset)
    L = get_attribute(filename, groupname, "L")
    newLayer = layer(x=prob, y=dataend, Geom.line, color = ["L = $L"])
end

f()