using HDF5
using Plots

function read_data(group)
    file = h5open(joinpath(@__DIR__, "supplementary_fig2_data.h5"), "r")
    y = read(file[group], "y")
    x = read(file[group], "x")
    close(file)
    return x, y
end

function plot_fig(fig)
    labels = ["random feedback measurement"]
    label = labels[1]
    x, y = read_data(label)
    pl = plot(x, y, label = label, xaxis = :log10)
    return pl
end

plot_fig("c")
