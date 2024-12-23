using JLD2
using Plots
using OrderedCollections
using PlotAndSave
using QuantumTimes
using LaTeXStrings

function T_fit(dir_name, start_cut)
    path = @__DIR__
    plotinfo = load(joinpath(dirname(path), dir_name, "data.jld2"), "plotinfo")
    lines = OrderedDict()

    for (key, line) in plotinfo.lines
        f_s = floor(Int, length(line.x) * start_cut)
        smaller_than = findfirst(x -> x < 0.01, line.y)
        f_e = isnothing(smaller_than) ? length(line.x) : smaller_than

        x = line.x[f_s:f_e]
        y = line.y[f_s:f_e]
        int, slp, fit = exp_fit(x, y)
        lines[line.tag] = fit, -1/slp
    end
    return lines
end

function TS_fit()
    T_fit("Ts", 0.01)
end
function T1_fit()
    T_fit("T1", 0.2)
end
function T2_fit()
    T_fit("T2", 0.1)
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
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
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
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
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
    ylabel = L"$\langle \sigma_x \rangle \textrm{envelope}$"
    )
    fit_lines = T2_fit()
    for line in values(plotinfo.lines)
        d = parse(Float64, line.tag[5:end])
        plot!(pl, line.x, line.y, label = "Disorder = " * string(d))
        f, T = fit_lines[line.tag]
        plot!(pl, line.x, f, label = "", lc = :grey, ls = :dash)
    end

    # normal decoherence
    # t = collect(values(plotinfo.lines))[1].x
    # deph_rate = parse(Float64, plotinfo.parameters[:deco][1:end - 3])
    # U = parse(Float64, plotinfo.parameters[:U][1:end - 3])
    # deph_rate = deph_rate / U
    # plot!(pl, t, x -> exp(-x * deph_rate * 0.5), label = "decoherence")

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
    t_time_xs = [parse(Float64, key[5:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    lines = T1_fit()
    t_time_x1 = [parse(Float64, key[5:end]) for key in keys(lines)]
    T1 = [val[2] for val in values(lines)]

    lines = T2_fit()
    t_time_x2 = [parse(Float64, key[5:end]) for key in keys(lines)]
    T2 = [val[2] for val in values(lines)]

    mega_cons = (10^6 * 2 * pi)
    data_path = joinpath(dirname(@__DIR__), "T1")
    J = get_param(data_path, :J) * mega_cons

    pl = plot(t_time_xs * J / mega_cons, Ts / J * 10^6,
    xlabel = L"$\textrm{Disorder}\ [\textrm{MHz}]$",
    ylabel = L"$\textrm{Survival\ Times}\ [\mu \textrm{s}]$",
    label = L"$T_s$",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9,
    yaxis = :log10
    )
    plot!(pl, t_time_x1 * J / mega_cons, T1 / J * 10^6, label = L"$T_1$")
    plot!(pl, t_time_x2 * J / mega_cons, T2 / J * 10^6, label = L"$T_2$")

    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.png"))
    savefig(pl, joinpath(@__DIR__, "combine_fit_plot.svg"))

    lines = []
    traj = get_traj(data_path)
    push!(lines, LineInfo(t_time_xs, Ts, traj, "Ts"))
    push!(lines, LineInfo(t_time_x1, T1, traj, "T1"))
    push!(lines, LineInfo(t_time_x2, T2, traj, "T2"))
    plotinfo = loaddata(data_path)
    savedata(@__DIR__, lines...; xlabel = "Disorder", ylabel = "Survival Times", copy_code = false, plotinfo.parameters...)
    return pl
end

function ts_proportional_check()
    lines = TS_fit()
    t_time_xs = [parse(Float64, key[5:end]) for key in keys(lines)]
    Ts = [val[2] for val in values(lines)]
    
    pl = plot(t_time_xs, Ts,
    xlabel = "J",
    ylabel = "Survival Times",
    label = "Ts",
    legend = :topleft,
    dpi = 300,
    titlefontsize = 9
    )
    L = 3
    U = 250
    f_ts(J) = 5/3*(L - 1) * 1 / (2 * (J/U)^2)
    plot!(pl, t_time_xs, f_ts, label = "5/3(L-1)U/(2J^2)")
    return pl
end

#Plots.scalefontsizes(1.5)
ts_plot()
t1_plot()
t2_plot()
combine_fit_plot()
println("Trajectories: $(get_traj(joinpath(dirname(@__DIR__), "T1")))")
# Plots.scalefontsizes();
# ts_proportional_check()