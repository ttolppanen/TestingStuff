using MIPTDataProcessing
using Gadfly

function getlayer(groupname)
    filename = "mipt_same_results"
    dt = get_attribute(filename, groupname, "dt")
    prob = get_probabilities(filename, groupname)
    # prob = filter(p -> p / dt > 0.49 && p / dt < 3.01, prob)
    prob = filter(p -> p / dt > 0.49 && p / dt < 2.01, prob)
    dataend = get_data_end_mean(filename, groupname, prob, "entanglement_half")
    L = get_attribute(filename, groupname, "L")
    return layer(x=prob, y=dataend, Geom.line, color = ["L = $L"])
end

function collapse_test()
    filename = "mipt_same_results"
    # data = ["mipt_19", "mipt_6", "mipt_20"]
    data = ["mipt_21", "mipt_11", "mipt_22"]
    for g in data
        println(get_attributes_string(filename, g; addnewline = false))
    end
    layers = [getlayer(g) for g in data]
    title = get_attributes_string(filename, data[1])
    p = plot(layers..., Guide.title(title), Guide.xlabel("p"), Guide.ylabel("S"), Coord.cartesian(xmin=0.04, xmax=0.16))
    p |> SVG("test.svg", 6inch, 4inch)
end

# f()
# accuracy_test()
collapse_test()