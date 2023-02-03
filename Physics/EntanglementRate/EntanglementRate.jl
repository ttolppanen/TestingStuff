using MIPTDataProcessing
using HDF5
using Gadfly

function f()
    filename = "mipt_same_results"
    groups = ["mipt_21", "mipt_11", "mipt_22"]
    for g in groups
        println(get_attributes_string(filename, g; addnewline = false))
    end
    layers = [getlayer(g, 0.2) for g in groups]
    title = get_attributes_string(filename, groups[1])
    p = plot(layers..., Guide.title(title), Guide.xlabel("t"), Guide.ylabel("(S(t + dt) - S(t)) / dt"), style(line_width=0.1mm))
    p |> SVG("test.svg", 6inch, 4inch)
end

function getlayer(groupname, p)
    filename = "mipt_same_results"
    data = get_data(filename, groupname, p, "entanglement_half")
    timesteps = size(data)[1]
    traj = size(data)[2]
    println(traj)
    y = zeros(timesteps - 1)
    for i in eachindex(y)
        y[i] = sum(data[i + 1, 1:traj] .- data[i, 1:traj]) / traj
    end
    dt = get_attribute(filename, groupname, "dt")
    L = get_attribute(filename, groupname, "L")
    x = [i * dt for i in 1:timesteps-1]
    y = y ./ dt
    newLayer = layer(x=x, y=y, Geom.smooth(smoothing = 0.01), color = ["p = $p, L = $L"])
end

f()