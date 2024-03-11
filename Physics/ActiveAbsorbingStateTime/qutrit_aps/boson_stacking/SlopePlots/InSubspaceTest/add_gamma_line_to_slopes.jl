using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections

function lin_reg_slope(x, y) # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    return (A \ y)[1]
end

function add_rabi(U, J, dt)
    w = sqrt(U^2 + 2*(J)^2) / (2 * pi)
    p = w * 2 * dt
    println(p)
    if U == 7
        color = :blue
        tag = "U = 7"
    elseif U == 9
        color = :green
        tag = "U = 9"
    elseif U == 10
        color = :gold
        tag = "U = 10"
    elseif U == 11
        color = :pink
        tag = "U = 11"
    else
        color = :cyan
        tag = "U = 13"
    end
    vline!([p], color = color, label = tag)
end

function combine_slope_plots()
    path = joinpath(@__DIR__, "L12")
    pl = plot(xlabel = "p", 
    ylabel = "Exp decay slopes", 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright
    )
    for U in [7, 9, 10, 11, 13]
        data_path = joinpath(path, "U$U")
        plotinfo = load(joinpath(data_path, "data.jld2"), "plotinfo")
        ps = []
        slopes = []
        for line in values(plotinfo.lines)
            p = parse(Float64, line.tag[5:end])
            gamma = p / plotinfo.parameters[:dt]
            x = line.x .* gamma
            y = log10.(line.y)
            line_start_index = Int(floor(length(x) / 4))
            line_end_index = Int(floor(length(x) * 3 / 4))
            slope = lin_reg_slope(x[line_start_index : line_end_index], y[line_start_index : line_end_index])
            push!(ps, p)
            push!(slopes, slope)
        end
        traj = collect(values(plotinfo.lines))[1].traj
        title!(plotinfo.title * ", traj = $traj")
        plot!(pl, ps, slopes, label = "U = $U")
        add_rabi(U, plotinfo.parameters[:J], plotinfo.parameters[:dt])
    end
    path = joinpath(path, "slopes with msr gamma.png")
    savefig(pl, path)
    return pl
end

combine_slope_plots()