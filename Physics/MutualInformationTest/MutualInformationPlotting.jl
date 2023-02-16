using MIPTDataProcessing
using Gadfly

function f()
    filename = "mipt_same_results"
    groups = ["mipt_19", "mipt_6", "mipt_20"]
    for g in groups
        println(get_attributes_string(filename, g; addnewline = false))
    end
    layers = [getlayer(g, "entanglement_half") for g in groups]
    title = get_attributes_string(filename, groups[1])
    p = plot(layers..., Guide.title(title), Guide.xlabel("Gamma"))
    p |> SVG("test.svg", 8inch, 6inch)
end

function getlayer(groupname, dataset)
    filename = "mipt_same_results"
    dt = get_attribute(filename, groupname, "dt")
    prob = get_probabilities_dataset(filename, groupname, dataset)
    prob = filter(p -> p > 0 && p < 0.05, prob)
    # prob = filter(p -> p / dt > 0.49 && p / dt < 2.01, prob)
    dataend = get_data_end_mean(filename, groupname, prob, dataset)
    L = get_attribute(filename, groupname, "L")
    newLayer = layer(x=prob ./ dt, y=dataend, Geom.line, color = ["L = $L"])
end

function collapse_test()
    filename = "mipt_same_results"
    dt02 = ["mipt_19", "mipt_6", "mipt_20"]
    dt08 = ["mipt_21", "mipt_11", "mipt_22"]
    for g in vcat(dt02, dt08)
        println(get_attributes_string(filename, g; addnewline = false))
    end
    layers = [getlayer(g) for g in vcat(dt02, dt08)]
    title = get_attributes_string(filename, dt02[1])
    p = plot(layers..., Guide.title(title * ", L = {6, 8, 10}"), Guide.xlabel("Gamma"))
    p |> SVG("test.svg", 6inch, 4inch)
end

f()
# collapse_test()