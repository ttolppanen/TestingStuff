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

function qubit_error_fix_U(sign, d, L, param, dw, dU, ranges, perm_mat)
    time = [0.0]

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
            if i != j && !((n_i == [1, 0] && n_j == [0, 1]) || (n_i == [0, 1] && n_j == [1, 0]))
                H[i, j] = 0
            end
        end
    end

    out(state, id) = begin
        time[1] += param[:dt]
        
        U = sparse(exp(sign * 1im * Matrix(H) * time[1]))
        U = kron(U, sparse(I, d^(param[:L] - L), d^(param[:L] - L)))
        
        return subspace_split(U, ranges, perm_mat)[id] * state[id], id
    end
    return out
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

function initial_parameters()
    param = OrderedDict()
    param[:gamma] = 0.0001
    param[:dt] = 0.02
    param[:t] = 90
    param[:d] = 3
    param[:L] = 2
    param[:state] = ValueInfo(bosonstack(1, param[:L], 1; d = param[:d]), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:w] = 1000
    param[:dw] = 0
    param[:U] = 50
    param[:J] = 1

    d = param[:d]; L = param[:L]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    traj = 1
    n_plots = 5
    # -----------

    d = param[:d]; L = param[:L]
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    n_op = [subspace_split(singlesite_n(d, L, i), ranges, perm_mat) for i in 1:L]
    obs = [(state, id) -> expval(state[id], n_op[i][id]) for i in 1:L]

    dt = param[:dt]; t = param[:t]
    pa_args, pa_out = pre_alloc_threads(zeros(length(obs) + 2, length(0:dt:t)), traj, deepcopy(state))

    out = []

    r_f(out, work_vector) = begin
        dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
        dU = 2 * dw .- (2 * param[:w] - param[:U])
        H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
        Uop = subspace_split(exp(-1im * dt * Matrix(H)), ranges, perm_mat)

        reverse_obs = qubit_error_fix_U(1, d, 2, param, dw, dU, ranges, perm_mat)
        reverse_back_obs = qubit_error_fix_U(-1, d, 2, param, dw, dU, ranges, perm_mat)

        first_obs = (state, id) -> begin
            new_state, id = reverse_obs(state, id)
            state[id] .= new_state
            return 100
        end
        last_obs = (state, id) -> begin
            new_state, id = reverse_back_obs(state, id)
            state[id] .= new_state
            return 100
        end

        effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
        exactevolve(state, work_vector, initial_id, Uop, dt, t, first_obs, obs... , last_obs; effect!, out)
    end

    plots = []
    for _ in 1:n_plots  
        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        map = heatmap(1:L, 0:dt:t, transpose(r[2:(end-1), :]), xlabel = "site", ylabel = "t", legend = :none)
        push!(plots, map)
    end

    title = PlotAndSave.parameters_to_title(; param...)
    pl = plot(plots..., layout = (1, n_plots), size = (1500, 1500), plot_title = title, plot_titlefont = 10)
    savefig(pl, joinpath(@__DIR__, "heatmap.png"))
end

@time f()