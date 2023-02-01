using MIPTDataProcessing
using Gadfly

function f()
    filename = "mipt_same_results"
    dt02 = "mipt_6"
    dt04 = "mipt_10"
    dt08 = "mipt_11"
    println(get_attributes_string(filename, dt02; addnewline = false))
    println(get_attributes_string(filename, dt04; addnewline = false))
    println(get_attributes_string(filename, dt08; addnewline = false))
    layers = [getlayer(dt02), getlayer(dt04), getlayer(dt08)]
    title = get_attributes_string(filename, dt02)
    p = plot(layers..., Guide.title(title), Guide.xlabel("Gamma"))
    p |> SVG("test.svg", 6inch, 4inch)
end

function getlayer(groupname)
    filename = "mipt_same_results"
    dt = get_attribute(filename, groupname, "dt")
    prob = get_probabilities(filename, groupname)
    # prob = filter(p -> p / dt > 0.49 && p / dt < 3.01, prob)
    prob = filter(p -> p / dt > 0.49 && p / dt < 2.01, prob)
    dataend = get_data_end_mean(filename, groupname, prob, :entanglement_half)
    newLayer = layer(x=prob ./ dt, y=dataend, Geom.line, color = ["dt = $dt"])
end

function accuracy_test()
    filename = "mipt_same_results"
    dt08 = "mipt_11"
    dt08k10 = "mipt_12"
    println(get_attributes_string(filename, dt08; addnewline = false))
    println(get_attributes_string(filename, dt08k10; addnewline = false))
    layers = [getlayer(dt08), getlayer(dt08k10)]
    title = get_attributes_string(filename, dt08k10)
    p = plot(layers..., Guide.title(title), Guide.xlabel("Gamma"))
    p |> SVG("test.svg", 6inch, 4inch)
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
    p = plot(layers..., Guide.title(title * ", L = {6, 8, 10}"), Guide.xlabel("p"))
    p |> SVG("test.svg", 6inch, 4inch)
end

# f()
# accuracy_test()
collapse_test()