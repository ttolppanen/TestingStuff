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
    dt, msr_rate, feedback::Vector; skip_subspaces = [])
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
        for site in [length(msr_op[current_subspace_id])]
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

function initial_parameters(L)
    param = OrderedDict()
    param[:gamma] = 27
    param[:dt] = 0.02
    param[:t] = 90
    param[:d] = 3
    param[:state] = ValueInfo(bosonstack(1, L, 1; d = param[:d]), "|100...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = ValueInfo(param[:w] .+ (rand(L) .- 0.5) * 10, "10")
    param[:U] = 50
    param[:dU] = ValueInfo(2 * param[:dw].val .- (2 * param[:w] - param[:U]), "10")
    param[:J] = 1

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
    traj = 1000
    Ls = [3,4,5,6,7]
    # -----------
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        initial_id = find_subspace(param[:state].val, indices)
        state = subspace_split(param[:state].val, ranges, perm_mat)

        msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
        feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

        op = param[:state].val * param[:state].val'
        op = subspace_split(op, ranges, perm_mat)
        obs = (state, id) -> expval(state[id], op[id])

        dt = param[:dt]; t = param[:t]
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, deepcopy(state))

        r_f(out, work_vector) = begin
            U = H_from_param(L, param, ranges, perm_mat)
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            exactevolve(state, work_vector, initial_id, U, dt, t, obs; effect!, out, save_only_last = true)
        end
        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        
        path = joinpath(@__DIR__, "$L")
        makehistogram(path, map(x -> x[1, 1], pa_out); xlabel = "initial state overlap", ylabel = "count", copy_code = false, bins = 10, param...)
    end
end

function H_from_param(L, param, ranges, perm_mat)
    dw = param[:w] .+ (rand(L) .- 0.5) * 10
    dU = 2 * dw .- (2 * param[:w] - param[:U])
    H = bosehubbard(param[:d], L; w = dw, U = dU, J = param[:J])
    U = subspace_split(exp(-1im * param[:dt] * Matrix(H)), ranges, perm_mat)
    return U
end

@time f()