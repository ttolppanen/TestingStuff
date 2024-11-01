using JLD2
using Plots
using OrderedCollections
using PlotAndSave

function T1(x, y) 
    # linreg for slope
    # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    res = (A \ log.(y))
    slope = res[1]
    c = res[2]
    return x -> exp.(slope * x .+ c), -1 / slope 
end

function T_fit(dir_name, start_cut; end_cut = 1.0)
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), dir_name, "data.jld2"), "plotinfo")
    lines = OrderedDict()

    for (key, line) in plotinfo.lines
        f_s = floor(Int, length(line.x) * start_cut)
        smaller_than = findfirst(x -> x < 0.05, line.y)
        f_e = isnothing(smaller_than) ? Int(floor(length(line.x) * end_cut)) : smaller_than

        x = line.x[f_s:f_e]
        y = line.y[f_s:f_e]
        lines[line.tag] = T1(x, y)
    end
    return lines
end

function find_TS_from_less_than(less_than = 0.01)
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "Ts", "data.jld2"), "plotinfo")
    Ts_times = []

    for (key, line) in plotinfo.lines
        smaller_than = findfirst(x -> x < less_than, line.y)
        Ts = isnothing(smaller_than) ? line.x[end] : line.x[smaller_than]
        push!(Ts_times, Ts)
    end
    return Ts_times
end

function TS_fit()
    T_fit("Ts", 0.01)
end
function T1_fit()
    T_fit("T1", 0.05)
end
function T2_fit()
    T_fit("T2", 0.05)
end

function ts_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "Ts", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "Array Leakage Population"
    )
    fit_lines = TS_fit()
    for line in values(plotinfo.lines)
        J = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "J = " * string(J))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    
    savefig(pl, joinpath(path, "Ts.png"))
    savefig(pl, joinpath(path, "Ts.svg"))
    return pl
end

function t1_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T1", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "First Site Population",
    legend = :bottomleft
    )
    fit_lines = T1_fit()
    for line in values(plotinfo.lines)
        J = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "J = " * string(J))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end
    plot!(pl, [0], [0], label = "fit", lc = :grey)
    savefig(pl, joinpath(path, "T1.png"))
    savefig(pl, joinpath(path, "T1.svg"))
    return pl
end

function t2_plot()
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), "T2", "data.jld2"), "plotinfo")
    pl = plot(dpi = 300,
    xlabel = "time", 
    ylabel = "<σₓ> envelope"
    )
    fit_lines = T2_fit()
    for line in values(plotinfo.lines)
        J = parse(Int64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "J = " * string(J))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end

    # normal decoherence
    t = plotinfo.lines["J = 5"].x
    deph_rate = parse(Float64, plotinfo.parameters[:deco][1:end - 3])
    U = parse(Float64, plotinfo.parameters[:U][1:end - 3])
    deph_rate = deph_rate / U
    plot!(pl, t, x -> exp(-x * deph_rate * 0.5), label = "decoherence")

    plot!(pl, [0], [0], label = "fit", lc = :grey)
    savefig(pl, joinpath(path, "T2.png"))
    savefig(pl, joinpath(path, "T2.svg"))
    return pl
end



function make_fit_line(dir, cut_size)
    path = joinpath(dirname(@__DIR__), dir)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")

    lines = values(plotinfo.lines)
    t_times = []
    t_time_x = []
    for line in lines
        t_x = parse(Int64, line.tag[5:end])
        fit_line_start = floor(Int, length(line.x) * cut_size)
        t_line, t_time = T1(line.x[fit_line_start:end], line.y[fit_line_start:end])

        push!(t_times, t_time)
        push!(t_time_x, t_x)
    end
    return t_time_x, t_times
end

function combine_fit_plot()
    lines = TS_fit()
    t_time_xs = [parse(Int64, key[5:end]) for key in keys(lines)]
    Ts = find_TS_from_less_than(0.1)
    
    lines = T1_fit()
    t_time_x1 = [parse(Int64, key[5:end]) for key in keys(lines)]
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    t_time_x2 = [parse(Int64, key[5:end]) for key in keys(lines)]
    T2 = [val[2] for val in values(lines)]

    pl = plot(t_time_xs, Ts,
    xlabel = "J",
    ylabel = "Survival Times",
    label = "Tₛ",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9
    )
    plot!(pl, t_time_x1, T1, label = "T₁")
    plot!(pl, t_time_x2, T2, label = "T₂")

    plotinfo = load(joinpath(dirname(@__DIR__), "T1", "data.jld2"), "plotinfo")
    to_flt(x) = parse(Float64, match(r"\d+", x).match)
    U = to_flt(plotinfo.parameters[:U])
    diss = to_flt(plotinfo.parameters[:diss]) * 10^(-3) / U
    hline!(pl, [1 / diss], ls = :dash, label = "inherit dissipation")
    deco = to_flt(plotinfo.parameters[:deco]) * 10^(-3) / U
    hline!(pl, [2 / deco], ls = :dash, label = "inherit decoherence")

    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))

    lines = []
    traj = get_traj(plotinfo)
    push!(lines, LineInfo(t_time_xs, Ts, traj, "Ts"))
    push!(lines, LineInfo(t_time_x1, T1, traj, "T1"))
    push!(lines, LineInfo(t_time_x2, T2, traj, "T2"))
    savedata(@__DIR__, lines...; xlabel = "J", ylabel = "Survival Times", copy_code = false, plotinfo.parameters...)
    return pl
end

#Plots.scalefontsizes(1.5)
ts_plot()
t1_plot()
t2_plot()
combine_fit_plot();
Plots.scalefontsizes();