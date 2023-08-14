using Distributed
using JLD2
using Gadfly
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function zeroone_orderparameter(L)
    d = 2 # not well defined for higher d
    order_parameter = singlesite_n(d, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        if isodd(i)
            order_parameter .+= singlesite_n(d, L, i)
        else
            order_parameter .+= singlesite([1 0; 0 1], L, 1) .- singlesite_n(d, L, i)
        end
    end
    order_parameter ./= L
    return state -> expval(state, order_parameter)
end

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = allone(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 0) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    # calc_order = zeroone_orderparameter(L)
    op = nall(d, L)
    calc_order(state) = expval(state, op) ./ L
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measuresitesrandomly!(state, msrop, p)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    dst_y = r[1, :]
    # r = solvetrajectories(r_f, traj; paral = :distributed)
    # for p in prob
    #     effect!(state) = measuresitesrandomly!(state, msrop, p; kwargs...)
    #     r_f() = krylovevolve(state0, H, dt, t, k; effect!)
    #     r = solvetrajectories(r_f, traj; paral = :distributed)
    #     out = traj_mean(r)
    #     push!(ent_y, out[1, 1])
    #     push!(fluc_y, out[2, 1])
    # end
    return [LineInfo(0:dt:t, dst_y, traj, "p = $p")]
end

function f(p)
    param = Dict(
        :dt => 0.02,
        :t => 10,
        :w => 0,
        :U => 0.14,
        :J => 1,
        :d => 3,
        :measurement => "std_n_feedback_00",
        :save_before_effect => false,
        :state => "ones",
        :alg => "krylov",
        :k => 6,
        :L => 8,
        :p => p,
        :traj => 100
    )
    dst_layers = krylov(param[:d], param[:L], param[:dt], param[:t], param[:traj], param[:p]; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
    println("Saving Data")
    path = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L$(param[:L])_p$(param[:p])")
    # addtrajectories(path, ent_layers...)
    makeplot(path, dst_layers...; xlabel = "t", ylabel = "N / L", param...)
end

function combine_figures()
    new_path = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_combined")
    path1 = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_p0.02")
    path2 = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_p0.04")
    combineplots(new_path, path1, path2)
    for p in 0.06:0.02:0.2
        path1 = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_p$p")
        combineplots(new_path, new_path, path1)
    end
end

function plot_to_log()
    path = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_combined")
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    lines = []
    for line in values(plotinfo.lines)
        log_line = LineInfo(log.(line.x), log.(line.y), line.traj, line.tag)
        push!(lines, log_line)
    end
    path = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_combined_log")
    makeplot(path, lines...; xlabel = "log(t)", ylabel = "log(N / L)", plotinfo.parameters...)
end

function change_log_plot_bounds()
    path = joinpath(@__DIR__, "order_parameter_test", "abs_d3_tFunc_L8_combined")
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    layers = []
    for line in values(plotinfo.lines)
        log_layer = layer(x=line.x, y=line.y, Geom.line, color = [line.tag])
        push!(layers, log_layer)
    end
    path = joinpath(path, "log plot.svg")
    p = plot(layers..., Coord.cartesian(ymin = -3), Scale.x_log, Scale.y_log, Guide.title(plotinfo.title), Guide.xlabel(plotinfo.xlabel), Guide.ylabel(plotinfo.ylabel))
    p |> SVG(path, 6inch, 4inch)
end

# for p in 0.02:0.02:0.2
#     @time f(p)
# end
combine_figures()
plot_to_log()
change_log_plot_bounds()