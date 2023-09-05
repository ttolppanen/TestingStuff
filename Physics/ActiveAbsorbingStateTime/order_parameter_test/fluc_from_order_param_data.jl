using JLD2
using PlotAndSave
using Statistics
using Plots

function make_figure()
    path_to_data = joinpath(@__DIR__, "tFunc_smaller_prob_L10")
    plotinfo = load(joinpath(path_to_data, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = "SSFluc_tuure",
    titlefontsize = 9, 
    dpi = 300)
    y_smooth_amount = 100
    for line in values(plotinfo.lines)
        x_ticks = length(line.x)
        y_ss = mean(line.y[Int(floor(x_ticks / 2)) : end])
        y = (line.y .- y_ss).^2
        y = [mean(y[i : i + y_smooth_amount]) for i in 1:(x_ticks - y_smooth_amount)]
        plot!(pl, line.x[1:end - y_smooth_amount], y, label = line.tag)
    end
    path = joinpath(path_to_data, "tuure_fluc.png")
    savefig(pl, path)
end

function make_figure_log()
    path_to_data = joinpath(@__DIR__, "tFunc_smaller_prob_L10")
    plotinfo = load(joinpath(path_to_data, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = "SSFluc_tuure",
    titlefontsize = 9, 
    dpi = 300,
    yaxis = :log10)
    y_smooth_amount = 100
    for line in values(plotinfo.lines)
        x_ticks = length(line.x)
        y_ss = mean(line.y[Int(floor(x_ticks / 2)) : end])
        y = (line.y .- y_ss).^2
        y = [mean(y[i : i + y_smooth_amount]) for i in 1:(x_ticks - y_smooth_amount)]
        plot!(pl, line.x[1:end - y_smooth_amount], y, label = line.tag)
    end
    path = joinpath(path_to_data, "tuure_fluc_log.png")
    savefig(pl, path)
end

make_figure()
make_figure_log()