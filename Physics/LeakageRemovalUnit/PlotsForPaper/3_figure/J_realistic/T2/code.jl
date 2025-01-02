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
    param[:T] = ValueInfo(0.1, "0.1K")
    param[:deco] = ValueInfo(2 / (10.0 * 10^(-6) * U), "10.0µs") # 2 due to qubit-n op equivalence
    param[:diss] = ValueInfo(1 / (16.7 * 10^(-6) * U), "16.7µs")
    param[:w] = ValueInfo(7500 * mega_const / U, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / U, "100MHz")
    param[:gamma] = 3/2 * (J / U)^2 # this is 3/2 J^2/U, in units of U

    param[:t] = 1000 * U / (5 * J) # J was 5Mhz
    param[:dt] = 0.1
    param[:d] = 3
    param[:L] = 3

    state = complex(spzeros(param[:d]))
    state[1] = 1.0
    state[2] = 1.0
    normalize!(state)
    param[:state] = ValueInfo(state, "|0>+|1> + thermal state")
    
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6

    param[:effect] = "rand start, msr_n last site -> fb -> 0"
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

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    fb_op = [1 1 1; 0 0 0; 0 0 0]
    feedback = [singlesite(fb_op, L, i) for i in 1:L]

    diss_op = [sqrt(param[:diss].val) * singlesite_a(d, L, i) for i in 1:L]
    deco_op = [sqrt(param[:deco].val) * singlesite_n(d, L, i) for i in 1:L]
    H_dd = complex(spzeros(d^L, d^L))
    for op in vcat(diss_op, deco_op)
        H_dd .+= -1.0im / 2 * op' * op
    end

    sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
    op = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, op)

    for J_i in Js
        param = initial_parameters(J_i)
        d = param[:d]; L = param[:L]
        dt = param[:dt]; t = param[:t]

        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - 1)
            H = bosehubbard(d, L; w = dw, U = dU, J = J_i / param[:U].val)
            H .+= H_dd

            t_T = param[:T].val
            t_w = param[:w].val * param[:U].val
            t_U = param[:U].val
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)

            effect! = periodic_msr_feedback!_function(msrop, param[:dt], param[:gamma], feedback, [param[:L]])
            U = exp(-1.0im * Matrix(H) * param[:dt])
            exactevolve(state, U, vcat(diss_op, deco_op), dt, t, obs; effect!)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "J = $(J_i / (10^6 * 2 * pi))MHz"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "sigma_x1 env", copy_code = false, add_traj = true,  param...)
    # addtrajectories(path, lines...)
end

@time f()