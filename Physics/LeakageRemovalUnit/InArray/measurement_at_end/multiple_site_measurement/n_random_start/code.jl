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

function measure_ends_effect!(msr_op::Vector, 
    dt, msr_rate, feedback::Vector, sites; skip_subspaces = [])
    time = [rand() / msr_rate]
    out(state, subspace_id) = begin
        time[1] += dt

        if time[1] < 1 / msr_rate
            return state, subspace_id::Int64
        end
        time[1] -= 1 / msr_rate

        if in(subspace_id, skip_subspaces)
            return state, subspace_id::Int64
        end
        current_subspace_id = subspace_id
        for site in sites
            msr_outcome = measuresite!(state[current_subspace_id], msr_op[current_subspace_id], site)
            new_id, op = feedback[current_subspace_id][site][msr_outcome]
            if new_id == current_subspace_id
                state[new_id] .= op * state[current_subspace_id]
            else
                mul!(state[new_id], op, state[current_subspace_id])
            end
            normalize!(state[new_id])
            current_subspace_id = new_id
        end
        return state, current_subspace_id::Int64
    end
    return out
end

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:t] = 90
    param[:d] = 3
    param[:L] = 7
    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")
    param[:alg] = "kv sbspc"
    param[:k] = 6
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = ValueInfo(param[:w] .+ (rand(param[:L]) .- 0.5) * 10, "10")
    param[:U] = 50
    param[:dU] = ValueInfo(2 * param[:dw].val .- (2 * param[:w] - param[:U]), "10")
    param[:J] = 1

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()

    # -----------
    traj = 2000
    gammas = vcat(0.01:0.01:0.1, 0.2:0.2:1, 1:30)
    site_targets = [[7], [6,7], [5,6,7], [4,5,6,7], [3,4,5,6,7]]
    # -----------

    d = param[:d]; L = param[:L];
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    anh = subspace_split(anhall(d, L), ranges, perm_mat)
    obs = (state, id) -> expval(state[id], anh[id])

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub(state, k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, pa_k)

    H = subspace_split(bosehubbard(d, L; w = param[:dw].val, U = param[:dU].val, J = param[:J]), ranges, perm_mat)
    
    lines = []
    for sites in site_targets
        out = []
        for gamma in gammas
            r_f(out, pa_k) = begin
                effect! = measure_ends_effect!(msrop, param[:dt], gamma, feedback, sites)
                krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs; effect!, out, save_only_last = true)
            end
            solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "$sites"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()
exit()