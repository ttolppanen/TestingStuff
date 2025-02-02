using Distributed
using JLD2
using Plots
using LinearAlgebra
using SparseArrays
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters(J)
    param = OrderedDict()
    mega_const = 10^6 * 2 * pi
    param[:U] = ValueInfo(250 * mega_const, "250MHz")
    U = param[:U].val
    param[:w] = ValueInfo(7500 * mega_const / U, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / U, "100MHz")
    param[:gamma] = 3/2 * (J / U)^2 # this is 3/2 J^2/U, in units of U

    param[:t] = 1000 * U / (5 * mega_const) # J was 5Mhz
    param[:dt] = 0.1
    param[:d] = 3
    param[:L] = 3

    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")
    
    param[:alg] = "exact dd"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]

    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(1)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Js = [15] .* (10^6 * 2 * pi)
    # -----------
    lines = []
    for J_i in Js
        param = initial_parameters(J_i)
        d = param[:d]; L = param[:L]
        
        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb
        
        op = anhall(d, L)
        obs = (state) -> expval(state, op)

        dt = param[:dt]; t = param[:t];

        dw = param[:w].val .+ ([1.0, 0.5, 0.0] .- 0.5) * param[:dw].val
        dU = 2 * dw .- (2 * param[:w].val - 1)
        H = bosehubbard(d, L; w = dw, U = dU, J = J_i / param[:U].val)

        # effect! = periodic_msr_feedback!_function(msrop, param[:dt], param[:gamma], feedback, [L])
        U = exp(-1.0im * Matrix(H) * param[:dt])
        r = exactevolve(param[:state].val, U, dt, t, obs)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "dt = $(param[:dt])"))
    end
    line = lines[1]
    plot(line.x, line.y, label = line.tag)
end

@time f()