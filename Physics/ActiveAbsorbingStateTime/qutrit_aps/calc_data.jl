using Distributed
using JLD2
using Plots
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

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = allone(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 0; 0 0 0; 0 0 1]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    n_d = nall(d, L) ./ L
    calc_order(state) = expval(state, n_d)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measurement_and_feedback!(state, msrop, p, fb_op)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_order; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    dst_y = r[1, :]
    return LineInfo(0:dt:t, dst_y, traj, "p = $p")
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
        :state => "allone",
        :alg => "krylov",
        :L => 8,
        :k => 6
    )
    traj = 200
    prob = 0.01:0.02:0.09
    dst_lines = []
    for p in prob
        order_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, p; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
        push!(dst_lines, order_line)
    end
    
    println("Saving Data")
    path = joinpath(@__DIR__, "boson_stacking" , "fb_00_10_22_U$(param[:U])_allone")
    addtrajectories(path, dst_lines...)
    # makeplot(path, dst_lines...; xlabel = "t", ylabel = "N / L", param...)
end

function linearize_plot(name)
    path = joinpath(@__DIR__, "boson_stacking", name)
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



# @time f()
# plot_to_log()
# linearize_plot("fb_00_10_22_U0.0")
# linearize_plot("fb_00_10_22_U1.0")
linearize_plot("fb_00_10_22_U10")
# linearize_plot("fb_00_10_22_U10_allone")
# linearize_plot_better()