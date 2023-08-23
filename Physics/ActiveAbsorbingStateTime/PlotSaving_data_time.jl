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
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    calc_order = zeroone_orderparameter(L)
    calc_ent(state) = entanglement(d, L, state, Int(floor(L/2)))
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measuresitesrandomly!(state, msrop, p)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order, calc_ent; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    dst_y = r[1, :]
    ent_y = r[2, :]
    # r = solvetrajectories(r_f, traj; paral = :distributed)
    # for p in prob
    #     effect!(state) = measuresitesrandomly!(state, msrop, p; kwargs...)
    #     r_f() = krylovevolve(state0, H, dt, t, k; effect!)
    #     r = solvetrajectories(r_f, traj; paral = :distributed)
    #     out = traj_mean(r)
    #     push!(ent_y, out[1, 1])
    #     push!(fluc_y, out[2, 1])
    # end
    return LineInfo(0:dt:t, dst_y, traj, "order, p = $p"), LineInfo(0:dt:t, ent_y, traj, "entanglement, p = $p")
end

function krylov_prob(d, L, dt, t, traj, prob; w, J, U, k)
    println("Generating")
    state0 = onezero(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(d, 0) : n_bosons_projector(d, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    # bra_zeroone = zeroone(d, L)'
    # calc_dst(state) = abs(bra_zeroone * state)
    calc_order = zeroone_orderparameter(L)
    H = bosehubbard(d, L; w, J, U)

    println("Solving trajectories")
    out = []
    r = solvetrajectories(r_f, traj; paral = :distributed)
    for p in prob
        effect!(state) = measuresitesrandomly!(state, msrop, p)
        r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!, save_only_last = true)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        push!(out, traj_mean(r)[1, 1])
    end
    return LineInfo(prob, out, traj, "L = $L")
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 10,
        :w => 0,
        :U => 0.14,
        :J => 1,
        :d => 2,
        :measurement => "std_n_feedback_01",
        :save_before_effect => false,
        :state => "onezero",
        :alg => "krylov",
        :L => 10,
        :entanglement => "half chain",
        :k => 6,
        :traj => 1000
    )
    prob = [0.01, 0.02, 0.04, 0.08]
    dst_lines = []
    for p in prob
        order_line, ent_line = krylov(param[:d], param[:L], param[:dt], param[:t], param[:traj], p; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line, ent_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "order_parameter_test", "ent_test")
    addtrajectories(path, dst_lines...)
    # makeplot(path, dst_lines...; xlabel = "t", param...)
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
    path = joinpath(@__DIR__, "order_parameter_test", "tFunc_new_L10")
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    layers = []
    for line in values(plotinfo.lines)
        log_layer = layer(x=line.x , y=line.y, Geom.line, color = [line.tag])
        push!(layers, log_layer)
    end
    path = joinpath(path, "log plot.svg")
    p = plot(layers..., Scale.y_log10, Guide.title(plotinfo.title), Guide.xlabel(plotinfo.xlabel), Guide.ylabel(plotinfo.ylabel))
    p |> SVG(path, 6inch, 4inch)
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

function linearize_plot()
    path = joinpath(@__DIR__, "order_parameter_test", "tFunc_new_L10")
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    layers = []
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        gamma = p / plotinfo.parameters[:dt]
        log_layer = layer(x=line.x .* gamma , y=line.y, Geom.line, color = [line.tag])
        push!(layers, log_layer)
    end
    path = joinpath(path, "linearized plot.svg")
    p = plot(layers..., Coord.cartesian(xmax = 5, ymin = -4), Scale.y_log10, Guide.title(plotinfo.title), Guide.xlabel("t*gamma"), Guide.ylabel(plotinfo.ylabel))
    p |> SVG(path, 6inch, 4inch)
end

@time f()
# plot_to_log()
# linearize_plot()