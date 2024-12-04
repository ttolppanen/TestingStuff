# this file had vanished, and was copied here from a previous commit
# titled "checked time evolution - 8ac916a9972601aaf0a32c769fa1c25e960d2b57"

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

function initial_parameters()
    param = OrderedDict()
    param[:dt] = 0.02
    param[:t] = 400
    param[:d] = 3
    param[:L] = 3
    param[:state] = ValueInfo(bosonstack(2, param[:L], 1), "|200...>")
    param[:alg] = "exact sbspc"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")

    d = param[:d]
    msrop = measurementoperators(nop(d), param[:L])
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, param[:L], i) for i in 1:param[:L]]

    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters()
    # -----------
    traj = 1000
    gammas = [0.03, 50]
    # -----------

    d = param[:d]; L = param[:L]
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    initial_id = find_subspace(param[:state].val, indices)
    state = subspace_split(param[:state].val, ranges, perm_mat)

    msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
    feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

    anh = subspace_split(anhall(d, L), ranges, perm_mat)
    obs = (state, id) -> expval(state[id], anh[id])

    dt = param[:dt]; t = param[:t]
    pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, deepcopy(state))
    
    lines = []
    for gamma in gammas
        r_f(out, work_vector) = begin
            dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
            dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
            H = bosehubbard(d, L; w = dw, U = dU, J = 1)
            Uop = subspace_split(exp(-1im * dt * Matrix(H)), ranges, perm_mat)
            effect! = periodic_sbspc_msr_feedback!_function(msrop, param[:dt], gamma, feedback, [L])
            exactevolve(state, work_vector, initial_id, Uop, dt, t, obs; effect!, out)
        end
        @time solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        y = r[1, :]
        push!(lines, LineInfo(0:dt:t, y, traj, "gamma = $gamma"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "anh_tot", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()