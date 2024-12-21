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
    param[:dt] = 0.02
    param[:gamma] = 0.03
    param[:t] = 1000
    param[:d] = 3
    param[:L] = L
    
    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")

    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
    mega_const = 10^6 * 2 * pi
    param[:J] = ValueInfo(5 * mega_const, "5MHz")
    param[:w] = ValueInfo(7500 * mega_const / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 * mega_const / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 * mega_const / param[:J].val, "250MHz")

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]

    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
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

        msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
        feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

        op = subspace_split(anhall(d, L), ranges, perm_mat)
        obs = (state, id) -> expval(state[id], op[id])

        dt = param[:dt]; t = param[:t]; k = param[:k]

        out = []
        r_f() = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)
            H = subspace_split(H, ranges, perm_mat)

            effect! = periodic_sbspc_msr_feedback!_function(msrop, param[:dt], param[:gamma], feedback, [L])
            krylovevolve(state, initial_id, H, dt, t, k, obs; effect!)
        end
        @time pa_out = solvetrajectories(r_f, traj; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "anh_tot", copy_code = false, add_traj = true, param...)
end

@time f()
