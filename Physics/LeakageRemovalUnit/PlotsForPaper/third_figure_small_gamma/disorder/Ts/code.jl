using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters()
    param = OrderedDict()
    param[:gamma] = 0.03
    param[:t] = 400
    param[:dt] = param[:t] / 1000
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"

    mega_const = 10^6 * 2 * pi
    param[:J] = ValueInfo(5 * mega_const, "5MHz")
    param[:w] = ValueInfo(7500 * mega_const / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 * mega_const / param[:J].val, "250MHz")

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]

    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    traj = 5000
    dws = range(1, 20, length = 100)
    # -----------
    lines = []
    d = param[:d]; L = param[:L]

    state = param[:state].val

    msrop = param[:effect].val.msr
    feedback = param[:effect].val.fb

    op = anhall(d, L)
    obs = (state) -> expval(state, op)

    dt = param[:dt]; t = param[:t]

    for dw_i in dws
        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)
            U = exp(-1.0im * Matrix(H) * param[:dt])
            effect! = periodic_msr_feedback!_function(msrop, param[:dt], param[:gamma], feedback, [param[:L]])
            exactevolve(state, U, dt, t, obs; effect!)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "anh_tot", copy_code = false, add_traj = true, param...)
    # addtrajectories(path, lines...)
end

@time f()