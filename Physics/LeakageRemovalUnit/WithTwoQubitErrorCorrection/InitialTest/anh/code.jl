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

function qubit_error_fix_U(d, L, param, dw, dU)
    H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
    n_i = zeros(L)
    n_j = zeros(L)
    for j in 1:d^L
        for site in 1:L
            n_i[site] = floor((j-1) / d^(site-1)) % d
        end
        for i in 1:d^L
            for site in 1:L
                n_j[site] = floor((i-1) / d^(site-1)) % d
            end
            if 2 in n_i || 2 in n_j
                H[i, j] = 0
            end
        end
    end
    U = sparse(exp(1im * Matrix(H) * param[:dt]))
    return kron(U, sparse(I, d^(param[:L] - L), d^(param[:L] - L)))
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
    param[:dt] = 0.02
    param[:t] = 90
    param[:d] = 3
    param[:L] = L
    param[:state] = ValueInfo(bosonstack(2, L, 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 10
    param[:U] = 50
    param[:J] = 1

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(4)
    # -----------
    traj = 2000
    Ls = [3,4,5,6,7]
    gammas = vcat(0.01:0.01:0.1, 0.2:0.2:1, 1:30)
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        initial_id = find_subspace(param[:state].val, indices)
        state = subspace_split(param[:state].val, ranges, perm_mat)

        msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
        feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

        anh = subspace_split(anhall(d, L), ranges, perm_mat)
        obs = (state, id) -> expval(state[id], anh[id])

        dt = param[:dt]; t = param[:t]
        pa_args, pa_out = pre_alloc_threads(zeros(1, 1), traj, deepcopy(state))

        dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
        dU = 2 * dw .- (2 * param[:w] - param[:U])
        H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
        Uop = subspace_split(exp(-1im * dt * Matrix(H)), ranges, perm_mat)

        error_U = qubit_error_fix_U(d, 2, param, dw, dU)
        error_U = subspace_split(error_U, ranges, perm_mat)
        fix_error(state, id) = begin
            state[id] .= error_U[id] * state[id]
            return state, id
        end

        out = []
        for gamma in gammas
            r_f(out, work_vector) = begin
                msr_effect! = measure_ends_effect!(msrop, param[:dt], gamma, feedback)
                effect!(state, id) = begin
                    state, id = fix_error(state, id)
                    return msr_effect!(state, id)
                end
                exactevolve(state, work_vector, initial_id, Uop, dt, t, obs; effect!, out, save_only_last = true)
            end
            solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
            r = traj_mean(pa_out)
            push!(out, r[1, 1])
        end
        push!(lines, LineInfo(gammas, out, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "gamma", ylabel = "anh", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()
exit()