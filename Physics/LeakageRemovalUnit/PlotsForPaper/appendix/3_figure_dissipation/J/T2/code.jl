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
using DSP
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
    param[:gamma] = 2 * (J / U)^2 # this is 2 J^2/U, in units of U

    param[:t] = 1000 * U / (5 * J) # J was 5Mhz
    param[:dt] = 0.1 * 20 * mega_const / J
    param[:d] = 3
    param[:L] = 3

    state = normalize(bosonstack(1, param[:L], 1; d = param[:d]) + allzero(param[:d], param[:L]))
    param[:state] = ValueInfo(state, "|10..> + |00..>")
    
    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
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
        
        sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
        op = singlesite(sigma_x, L, 1)
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
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "J = $(J_i / (10^6 * 2 * pi))MHz"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "J", ylabel = "sigma_x1 env", copy_code = false, add_traj = true,  param...)
    # addtrajectories(path, lines...)
end

@time f()