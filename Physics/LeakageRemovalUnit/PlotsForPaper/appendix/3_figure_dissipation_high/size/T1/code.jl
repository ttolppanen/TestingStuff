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
    
    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = param[:d]), "|100...>")

    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
    mega_const = 10^6 * 2 * pi
    param[:J] = ValueInfo(5 * mega_const, "5MHz")
    param[:w] = ValueInfo(7500 * mega_const / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 * mega_const / param[:J].val, "250MHz")
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
        
        initial_id = find_subspace(param[:state].val, indices)
        state = subspace_split(param[:state].val, ranges, perm_mat)

        op = subspace_split(singlesite_n(d, L, 1), ranges, perm_mat)
        obs = (state, id) -> expval(state[id], op[id])

        dt = param[:dt]; t = param[:t]; k = param[:k]

        diss_op = sqrt(param[:gamma]) * singlesite_a(d, L, L)
        H_dd = complex(spzeros(d^L, d^L))
        H_dd .+= -1.0im / 2 * diss_op' * diss_op

        a_relations = [-1 for _ in eachindex(ranges)]
        a_relations[1] = 0
        diss_op = subspace_split(diss_op, ranges, perm_mat, a_relations)

        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)
            H .+= H_dd
            H = subspace_split(H, ranges, perm_mat)

            krylovevolve(state, initial_id, H, [diss_op], [], dt, t, k, obs)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "n1", copy_code = false, add_traj = true, param...)
end

@time f()
