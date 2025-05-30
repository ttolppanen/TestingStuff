using HDF5
using Plots

function read_data(plot, group)
    file = h5open(joinpath(@__DIR__, "fig3_data.h5"), "r")
    y = read(file[plot][group], "y")
    x = read(file[plot][group], "x")
    close(file)
    return x, y
end

function plot_fig(fig)
    labels = ["ts", "t1", "t2"]
    label = labels[1]
    x, y = read_data(fig, label)
    pl = plot(x, y, label = label, yaxis = :log10)
    for label in labels[2:3]
        x, y = read_data(fig, label)
        plot!(pl, x, y, label = label)
    end
    return pl
end

plot_fig("c")
