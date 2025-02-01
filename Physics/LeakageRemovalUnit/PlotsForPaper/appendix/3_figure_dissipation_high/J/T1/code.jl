using Distributed
using JLD2
using Plots
using LinearAlgebra
using SparseArrays
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections
using SparseArrays

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
    param[:t] = 1000 * U / (5 * J) # J was 5Mhz
    param[:gamma] = 2
    param[:dt] = 1 / (10 * param[:gamma])
    param[:d] = 3
    param[:L] = 3

    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = param[:d]), "|100...>")

    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"
    return param
end

function f()
    param = initial_parameters(1)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Js = range(1, 20, length = 100) .* (10^6 * 2 * pi)
    # -----------
    lines = []
    for J_i in Js
        param = initial_parameters(J_i)
        d = param[:d]; L = param[:L]
        
        op = singlesite_n(d, L, 1)
        obs = (state) -> expval(state, op)

        dt = param[:dt]; t = param[:t];

        diss_op = sqrt(param[:gamma]) * singlesite_a(d, L, L)
        H_dd = complex(spzeros(d^L, d^L))
        H_dd .+= -1.0im / 2 * diss_op' * diss_op

        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - 1)
            H = bosehubbard(d, L; w = dw, U = dU, J = J_i / param[:U].val)
            H .+= H_dd
            U = exp(-1.0im * Matrix(H) * param[:dt])
            exactevolve(param[:state].val, U, [diss_op], dt, t, obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "J = $(J_i / (10^6 * 2 * pi))MHz"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "J", ylabel = "n_1", copy_code = false, add_traj = true, param...)
    # addtrajectories(path, lines...)
end

@time f()