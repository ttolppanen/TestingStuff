module LRUFigureColors

using Plots
using PlotAndSave
using OrderedCollections
using QuantumTimes

export color_small_msr
export color_high_msr
export color_periodic
export color_dissipation_removal
export color_t2
export color_t1
export color_ts
export color_ts_short_fig2
export color_ts_short

export T_fit
export t_small_scale

color_small_msr = colorant"#3c93c2" # blue
color_high_msr = colorant"#9ec9e2" # light blue
color_periodic = colorant"#00cd6c" # green
color_dissipation_removal = colorant"#9ccea7" # light green
color_t2 = colorant"#af58ba" # purple
color_t1 = colorant"#ffc61e" # yellow
color_ts = color_periodic
color_ts_short_fig2 = :lightgrey
color_ts_short = colorant"#06592A" # dark green

function T_fit(path, start_cut, end_cut = 1.0, y_factor = 1)
    plotinfo = loaddata(path)
    lines = OrderedDict()

    for (key, line) in plotinfo.lines
        f_s = floor(Int, length(line.x) * start_cut)
        smaller_than = findfirst(x -> x < 0.01, line.y)
        f_e = isnothing(smaller_than) ? floor(Int, length(line.x) * end_cut) : smaller_than

        x = line.x[f_s:f_e]
        y = line.y[f_s:f_e] * y_factor
        int, slp, fit = exp_fit(x, y)
        lines[line.tag] = fit, -1/slp
    end
    return lines
end

function t_small_scale(J, U)
    J = 2 * J^2 / U # 2 * J^2 / U
    return π / (2*J) # π / (2*J)
end
function t_small_scale(J, U, L)
    println("small scale time does not use L anymore")
    return t_small_scale(J, U)
end


# color_small_msr = RGB(245 /255, 76 /255, 181 /255)
# color_high_msr = RGB(91 /255, 245/255, 76/255)
# color_periodic = :black
# color_dissipation_removal = :black
# color_t2 = RGB(224 /255, 154 /255, 22 /255)
# color_t1 = RGB(224 /255, 203 /255, 22 /255)
# color_ts = color_periodic

# color_small_msr = :green
# color_high_msr = :red
# color_periodic = :blue
# color_dissipation_removal = :yellow
# color_t2 = :orange
# color_t1 = :yellow
# color_ts = color_periodic

end # module FigureColors
