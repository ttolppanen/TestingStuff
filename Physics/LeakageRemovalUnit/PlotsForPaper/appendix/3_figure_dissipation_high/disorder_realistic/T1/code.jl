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

function initial_parameters()
    param = OrderedDict()
    param[:t] = 400
    param[:d] = 3
    param[:L] = 3

    state = complex(spzeros(param[:d]))
    state[2] = 1.0
    param[:state] = ValueInfo(state, "|1> + thermal state")
    
    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"

    mega_const = 10^6 * 2 * pi
    param[:J] = ValueInfo(5 * mega_const, "5MHz")
    param[:w] = ValueInfo(7500 * mega_const / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 * mega_const / param[:J].val, "250MHz")
    param[:T] = ValueInfo(0.1, "0.1K")
    param[:diss] = ValueInfo(1 / (16.7 * 10^(-6) * param[:J].val), "16.7µs")
    param[:deco] = ValueInfo(2 / (10.0 * 10^(-6) * param[:J].val), "10.0µs") # 2 due to qubit-n op equivalence
    param[:gamma] = 2 * param[:U].val
    param[:dt] = 1 / (10 * param[:gamma])
    return param
end

function f()
    param = initial_parameters()
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    dws = range(1, 20, length = 100)
    # -----------
    lines = []

    d = param[:d]; L = param[:L]
    state = param[:state].val

    op = singlesite_n(d, L, 1)
    obs = (state) -> expval(state, op)

    dt = param[:dt]; t = param[:t];

    diss_op = [sqrt(param[:diss].val) * singlesite_a(d, L, i) for i in 1:(L - 1)]
    push!(diss_op, sqrt(param[:gamma]) * singlesite_a(d, L, L))
    deco_op = [sqrt(param[:deco].val) * singlesite_n(d, L, i) for i in 1:(L - 1)]
    H_dd = complex(spzeros(d^L, d^L))
    for op in vcat(diss_op, deco_op)
        H_dd .+= -1.0im / 2 * op' * op
    end
    
    for dw_i in dws
        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)
            H .+= H_dd
            U = exp(-1.0im * Matrix(H) * param[:dt])

            t_T = param[:T].val
            t_w = param[:w].val * param[:J].val
            t_U = param[:U].val * param[:J].val
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)
            exactevolve(state, U, vcat(diss_op, deco_op), dt, t, obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, add_traj = true, param...)
    # addtrajectories(path, line)
end

@time f()