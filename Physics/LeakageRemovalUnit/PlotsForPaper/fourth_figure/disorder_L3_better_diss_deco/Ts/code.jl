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
    param[:dt] = 0.02
    param[:gamma] = 0.03
    param[:t] = 400
    param[:d] = 3
    param[:L] = 3

    state = complex(spzeros(param[:d]))
    state[param[:d]] = 1.0
    param[:state] = ValueInfo(state, "|2> + thermal state")
    param[:T] = 0.1
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34
    J = 5 * 10^6 * 2 * pi # added 2 * pi 
    param[:J_value] = ValueInfo(J, "5MHz * 2pi")
    param[:T_scale] = ValueInfo(kb / (hbar * J), "kb/hbarJ")
    param[:diss] = 0.000106 # 300 mikrosecond
    param[:deco] = 0.000318 # 100 mikrosecond
    
    param[:alg] = "krylov sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:k] = 6
    param[:w] = 1500 # changed this from 1000
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
    traj = 1000
    dws = [5, 10, 15, 20]
    # -----------
    lines = []
    d = param[:d]; L = param[:L]
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    op = subspace_split(anhall(d, L), ranges, perm_mat)
    obs = (state, id) -> expval(state[id], op[id])

    dt = param[:dt]; t = param[:t]; k = param[:k]
    pa_k = PA_krylov_sub_dd(subspace_split(allzero(d, L), ranges, perm_mat), k)
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, pa_k)

    diss_op = [sqrt(param[:diss]) * singlesite_a(d, L, i) for i in 1:L]
    deco_op = [sqrt(param[:deco]) * singlesite_n(d, L, i) for i in 1:L]
    H_dd = complex(spzeros(d^L, d^L))
    for op in vcat(diss_op, deco_op)
        H_dd .+= -1.0im / 2 * op' * op
    end
    a_relations = [-1 for _ in eachindex(ranges)]
    a_relations[1] = 0
    diss_op = [subspace_split(op, ranges, perm_mat, a_relations) for op in diss_op]
    deco_op = [subspace_split(op, ranges, perm_mat) for op in deco_op]

    for dw_i in dws
        out = []
        r_f(out, pa_k) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * dw_i
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])

            t_T = param[:T]
            t_w = param[:w] * param[:J_value].val
            t_U = param[:U] * param[:J_value].val
            ts_sample = sample_transmon_thermal_state(d, L-1, t_T, t_w, t_U)
            state = kron(param[:state].val, ts_sample)

            initial_id = find_subspace(state, indices)
            state = subspace_split(state, ranges, perm_mat)

            H .+= H_dd
            H = subspace_split(H, ranges, perm_mat)
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            krylovevolve(state, initial_id, H, diss_op, deco_op, dt, t, k, pa_k, obs; effect!, out)
        end
        r_f(zeros(1, length(0:dt:t)), pa_k) # this is necessary for some reason?
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "dw = $dw_i"))
    end
    path = joinpath(@__DIR__)
    # makeplot(path, lines...; xlabel = "t", ylabel = "anh_tot", copy_code = false, param...)
    addtrajectories(path, lines...)
end

@time f()