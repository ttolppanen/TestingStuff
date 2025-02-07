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

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function initial_parameters(L)
    param = OrderedDict()
    param[:t] = 400
    param[:d] = 3
    param[:L] = L
    
    state = complex(spzeros(param[:d]))
    state[2] = 1.0
    param[:state] = ValueInfo(state, "|1> + thermal state")

    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
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
    param = initial_parameters(3)
    # -----------
    length(ARGS) > 0 ? traj = parse(Int, ARGS[1]) : traj = 100
    println("Number of trajectories: $traj")
    Ls = 2:6
    # -----------
    lines = []
    
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]; L = param[:L]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)

        op = subspace_split(singlesite_n(d, L, 1), ranges, perm_mat)
        obs = (state, id) -> expval(state[id], op[id])

        dt = param[:dt]; t = param[:t]; k = param[:k]
    
        diss_op = [sqrt(param[:diss].val) * singlesite_a(d, L, i) for i in 1:(L-1)]
        push!(diss_op, sqrt(param[:gamma]) * singlesite_a(d, L, L))
        deco_op = [sqrt(param[:deco].val) * singlesite_n(d, L, i) for i in 1:(L-1)]
        H_dd = complex(spzeros(d^L, d^L))
        for op in vcat(diss_op, deco_op)
            H_dd .+= -1.0im / 2 * op' * op
        end
        a_relations = [-1 for _ in eachindex(ranges)]
        a_relations[1] = 0
        diss_op = [subspace_split(op, ranges, perm_mat, a_relations) for op in diss_op]
        deco_op = [subspace_split(op, ranges, perm_mat) for op in deco_op]

        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)

            t_T = param[:T].val
            t_w = param[:w].val * param[:J].val
            t_U = param[:U].val * param[:J].val
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)

            initial_id = find_subspace(state, indices)
            state = subspace_split(state, ranges, perm_mat)

            H .+= H_dd
            H = subspace_split(H, ranges, perm_mat)
            krylovevolve(state, initial_id, H, diss_op, deco_op, dt, t, k, obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, add_traj = true, cut_data = true, param...)
end

@time f()
