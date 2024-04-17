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
    H = sparse(complex(zeros(9, 9)))
    H[d + 1, 2] = param[:J]
    H[2, d + 1] = param[:J]

    U = sparse(exp(1im * Matrix(H) * param[:dt]))
    sub_d = d^(param[:L]- L)
    return kron(U, sparse(I, sub_d, sub_d))
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
    param[:gamma] = 27
    param[:dt] = 0.5
    param[:t] = 400
    param[:d] = 3
    param[:L] = 5
    N = 1
    param[:state] = ValueInfo(bosonstack(N, param[:L], 1; d = param[:d]), "|$(N)00...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:extra_U] = "two qubit sigma_x"
    param[:w] = 1000
    param[:dw] = 1
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
    pa_args, pa_out = pre_alloc_threads(zeros(length(obs), length(0:dt:t)), traj, deepcopy(state))

    out = []

    r_f(out, work_vector) = begin
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

        msr_effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
        effect!(state, id) = begin
            state, id = fix_error(state, id)
            return msr_effect!(state, id)
        end
        exactevolve(state, work_vector, initial_id, Uop, dt, t, obs...; effect!, out)
    end

    plots = []
    for _ in 1:n_plots  
        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        map = heatmap(1:L, 0:dt:t, transpose(r), xlabel = "site", ylabel = "t", legend = :none)
        push!(plots, map)
    end

    title = PlotAndSave.parameters_to_title(; param...)
    pl = plot(plots..., layout = (1, n_plots), size = (1500, 1500), plot_title = title, plot_titlefont = 10)
    savefig(pl, joinpath(@__DIR__, "heatmap.png"))
end

@time f()