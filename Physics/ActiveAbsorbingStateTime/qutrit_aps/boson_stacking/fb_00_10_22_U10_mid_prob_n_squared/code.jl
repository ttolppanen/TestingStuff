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
    state0 = zeroone(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 0; 0 0 0; 0 0 1]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    calc_order = n_squared_orderparameter(d, L)
    calc_anharmonicity = anharmonicity_orderparameter(d, L)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measurement_and_feedback!(state, msrop, p, fb_op)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order, calc_anharmonicity; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    return LineInfo(0:dt:t, r[1, :], traj, "p = $p"), LineInfo(0:dt:t, r[2, :], traj, "p = $p")
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 10,
        :J => 1,
        :d => 3,
        :measurement => "std_n fb: 0 -> |0>, 1 -> |0>, 2 -> |2>",
        :save_before_effect => false,
        :state => "zeroone",
        :alg => "krylov",
        :L => 8,
        :k => 6
    )
    traj = 100
    prob = 0.01:0.025:0.21
    dst_lines = []
    anh_lines = []
    for p in prob
        order_line, anharmonicity_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, p; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line)
        push!(anh_lines, anharmonicity_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "boson_stacking" , "fb_00_10_22_U$(param[:U])_mid_prob_n_squared")
    # addtrajectories(path, dst_lines...)
    makeplot(path, dst_lines...; xlabel = "t", ylabel = "sum(N²ᵢ) / 2L", param...)
    linearize_plot(path)
    log_plot(path)

    path = joinpath(@__DIR__, "boson_stacking" , "fb_00_10_22_U$(param[:U])_mid_prob_anharmonicity")
    # addtrajectories(path, dst_lines...)
    makeplot(path, anh_lines...; xlabel = "t", ylabel = "sum(N(N - 1)) / L", param...)
    linearize_plot(path)
    log_plot(path)
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
    ylims = (10^-4, 1)
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



@time f()
