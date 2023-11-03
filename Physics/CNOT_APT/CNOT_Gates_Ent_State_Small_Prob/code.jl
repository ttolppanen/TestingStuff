using Distributed
using JLD2
using Plots
using LinearAlgebra
using SparseArrays
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

function anharmonicity_orderparameter(d, L)
    local_op = nop(d) * (nop(d) - I)
    out = singlesite(local_op, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        out .+= singlesite(local_op, L , i)
    end
    out ./= L
    return state -> expval(state, out)
end

function cnot_unitary(d, L)
    cnot_gate = sparse([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0])
    ⊗(a, b) = kron(a, b)
    single_site_gate(target) = sparse(I, d^(target-1), d^(target-1)) ⊗ cnot_gate ⊗ sparse(I, d^(L-target-1), d^(L-target-1))
    out = single_site_gate(L-1)
    for i in (L-2):-1:1
        out = single_site_gate(i) * out 
    end
    return out
end

function plus_minus_state(L) # not defined for d != 2
    plus = sparse([1/sqrt(2), 1/sqrt(2)])
    minus = sparse([1/sqrt(2), -1/sqrt(2)])
    sites = [isodd(i) ? plus : minus for i in 1:L]
    return complex(kron(sites...))
end

function krylov(d, L, dt, t, traj, p)
    println("Generating")
    state0 = plus_minus_state(L)
    msrop = measurementoperators(nop(d), L)
    feedback = [singlesite(n_bosons_projector(d, 0), L, i) for i in 1:L]
    n = nall(d, L) ./ L
    calc_order(s) = expval(s, n)
    H = cnot_unitary(d, L)

    effect!(state) = measurement_and_feedback!(state, msrop, p, feedback)
    r_f() = exactevolve(state0, H, dt, t, calc_order; effect!)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    return LineInfo(0:dt:t, r[1, :], traj, "p = $p")
end

function f()
    param = Dict(
        :dt => 0.1,
        :t => 30,
        :d => 2,
        :measurement => "std_n->0",
        :save_before_effect => false,
        :state => "|+>|->...", #there was a typo here, it said onezero
        :alg => "exact",
        :U => "cnot chain: C_12C_23...C_n,n+1",
        :L => 8
    )
    traj = 1000
    prob = 0.01:0.01:0.1
    dst_lines = []
    for p in prob
        order_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, p)
        push!(dst_lines, order_line)
    end
    
    println("Saving Data")
    path = @__DIR__
    # addtrajectories(path, dst_lines...)
    makeplot(path, dst_lines...; copy_code = false, xlabel = "t", ylabel = "N / L", param...)
    linearize_plot(path)
    log_plot(path)
    collapse_plot(path)
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
    ylims = (10^-5, 1)
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

@time f()
# collapse_plot(joinpath(@__DIR__, "high_prob"))