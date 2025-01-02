using Distributed
using JLD2
using Plots
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using SparseArrays

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters(L)
    param = OrderedDict()
    param[:dt] = 0.1
    param[:t] = 200
    param[:d] = 3

    state = complex(spzeros(param[:d]))
    state[param[:d]] = 1.0
    param[:state] = ValueInfo(state, "|2> + thermal state")

    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    mega_const = 10^6 * 2 * pi
    param[:J] = ValueInfo(5 * mega_const, "5MHz")
    param[:w] = ValueInfo(7500 * mega_const / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 * mega_const / param[:J].val, "250MHz")
    param[:T] = ValueInfo(0.1, "0.1K")
    param[:diss] = ValueInfo(1 / (16.7 * 10^(-6) * param[:J].val), "16.7µs")
    param[:deco] = ValueInfo(2 / (10.0 * 10^(-6) * param[:J].val), "10.0µs") # 2 due to qubit-n op equivalence


    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(2)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Ls = [3]
    gammas = exp.(x_ticks(log(10^-5), log(200), 200, 
        (log(0.009), log(0.04), 4), # more points between 0.009 and 0.04
        (log(3), log(10), 8))
        )
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        anh = anhall(d, L)
        obs = (state) -> expval(state, anh)

        dt = param[:dt]; t = param[:t];

        diss_op = [sqrt(param[:diss].val) * singlesite_a(d, L, i) for i in 1:L]
        deco_op = [sqrt(param[:deco].val) * singlesite_n(d, L, i) for i in 1:L]
        dd_op = vcat(diss_op, deco_op)
        H_dd = complex(spzeros(d^L, d^L))
        for op in dd_op
            H_dd .+= -1.0im / 2 * op' * op
        end

        out = []
        for gamma in gammas
            dt = min(param[:dt], 1 / (gamma))
            r_f() = begin
                dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
                dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
                H = bosehubbard(d, L; w = dw, U = dU, J = 1)
                H .+= H_dd
                Uop = exp(-1im * dt * Matrix(H))
                effect! = periodic_msr_feedback!_function(msrop, dt, gamma, feedback, [L])

                t_T = param[:T].val
                t_w = param[:w].val * param[:J].val
                t_U = param[:U].val * param[:J].val
                ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
                state = kron(param[:state].val, ts_sample)

                exactevolve(state, Uop, dd_op, dt, t, obs; effect!, save_only_last = true)
            end
            @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, add_traj = true, param...)
    # addtrajectories(path, lines...)
    xlog_plot(@__DIR__)
end

function xlog_plot(path)
    plotinfo = load(joinpath(path, "data.jld2"), "plotinfo")
    pl = plot(xaxis = :log10,
    title = plotinfo.title, 
    xlabel = plotinfo.xlabel, 
    ylabel = plotinfo.ylabel, 
    dpi = 300,
    titlefontsize = 7,
    legend = :outertopright
    # ylims = (10^-5, 1)
    )
    for line in values(plotinfo.lines)
        plot!(pl, line.x, line.y, label = line.tag)
    end
    path = joinpath(path, "log plot.png")
    savefig(pl, path)
    return pl
end

@time f()
# xlog_plot(@__DIR__)
