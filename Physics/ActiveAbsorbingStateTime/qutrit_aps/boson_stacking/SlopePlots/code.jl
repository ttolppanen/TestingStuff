using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

@everywhere function measurement_and_feedback!(state::AbstractVector{<:Number}, msr_op, msr_prob::Real, fb_op)
    for i in 1:length(msr_op)
        if rand(Float64) < msr_prob
            measuresite!(state, msr_op, i)
            state .= fb_op[i] * state
            normalize!(state)
        end
    end
end

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

function n_squared_orderparameter(d, L)
    n2 = nop(d).^2
    order_parameter = singlesite(n2, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        order_parameter .+= singlesite(n2, L, i)
    end
    order_parameter ./= (2 * L)
    return state -> expval(state, order_parameter)
end

function anharmonicity_orderparameter(d, L)
    local_op = nop(d) * (nop(d) - I)
    out = singlesite(local_op, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        out .+= singlesite(local_op, L , i)
    end
    out ./= L
    return state -> expval(state, out)
end

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = productstate(d, [2, 0, 2, 0, 2, 0, 2, 0])
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 0; 0 0 0; 0 0 1]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    n = nall(d, L) ./ L
    calc_order(s) = expval(s, n)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measurement_and_feedback!(state, msrop, p, fb_op)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    return LineInfo(0:dt:t, r[1, :], traj, "p = $p")
end

function f(U)
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => U,
        :J => 1,
        :d => 3,
        :measurement => "std_n fb: 0 -> |0>, 1 -> |0>, 2 -> |2>",
        :save_before_effect => false,
        :state => "|20202020>",
        :alg => "krylov",
        :L => 8,
        :k => 6
    )
    traj = 200
    prob = 0.01:0.01:0.1
    dst_lines = []
    for p in prob
        order_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, p; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "U$(param[:U])")
    addtrajectories(path, dst_lines...)
    # makeplot(path, dst_lines...; xlabel = "t", ylabel = "N / L", param...)
    linearize_plot(path)
    slope_plot(path)
end

function linearize_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(yaxis = :log10, 
    title = plotinfo.title, 
    xlabel = "t*gamma", 
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    xlims = (0, 30),
    ylims = (10^-1, 1)
    )
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        gamma = p / plotinfo.parameters[:dt]
        plot!(pl, line.x .* gamma, line.y, label = line.tag)
    end
    path = joinpath(path, "linearized plot.png")
    savefig(pl, path)
    return pl
end

function log_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(yaxis = :log10, 
    title = plotinfo.title, 
    xlabel = "t", 
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    ylims = (10^-5, 1)
    )
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y, label = line.tag)
    end
    path = joinpath(path, "log plot.png")
    savefig(pl, path)
    return pl
end

function test_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(yaxis = :log10, 
    title = plotinfo.title, 
    xlabel = "t/gamma",
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    xlims = (0, 30),
    ylims = (10^-4, 1)
    )
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        gamma = p / plotinfo.parameters[:dt]
        plot!(pl, line.x ./ gamma, line.y, label = line.tag)
    end
    path = joinpath(path, "linearized plot per gamma.png")
    savefig(pl, path)
    return pl
end

function lin_reg_slope(x, y) # https://jahoo.github.io/2021/01/11/simplest_linear_regression_example.html
    A = hcat(x, ones(length(x)))
    return (A \ y)[1]
end

function slope_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = "p", 
    ylabel = "Exp decay slopes", 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright
    )
    ps = []
    slopes = []
    for line in values(plotinfo.lines)
        p = parse(Float64, line.tag[5:end])
        gamma = p / plotinfo.parameters[:dt]
        x = line.x .* gamma
        y = log10.(line.y)
        slope = lin_reg_slope(x[1: Int(floor(end / 2))], y[1: Int(floor(end / 2))])
        push!(ps, p)
        push!(slopes, slope)
    end
    plot!(pl, ps, slopes)
    path = joinpath(path, "slope plot.png")
    savefig(pl, path)
    return pl
end

function collapse_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(title = plotinfo.title, 
    xlabel = "t", 
    ylabel = "t^δ * N / L", 
    dpi = 300,
    titlefontsize = 9,
    legend = :outertopright,
    xlims = (0, 30)
    )
    delta = 0.159
    dt = plotinfo.parameters[:dt]
    for line in values(plotinfo.lines)
        plot!(pl, line.x, [line.y[i] * (dt * (i-1))^delta for i in eachindex(line.y)], label = line.tag)
    end
    path = joinpath(path, "collapse plot.png")
    savefig(pl, path)
    return pl
end

@time f(7)
@time f(9)
@time f(10)
@time f(11)
@time f(13)
