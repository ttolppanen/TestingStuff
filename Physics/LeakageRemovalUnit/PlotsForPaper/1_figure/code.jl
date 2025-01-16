using Distributed
using JLD2
using Plots
using Plots.PlotMeasures
using LinearAlgebra
using QuantumStates
using QuantumOperators
using PlotAndSave
using QuantumTimeEvolution
using OrderedCollections
using LaTeXStrings

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function measure_ends_effect!(msr_op::Vector, 
    dt, msr_rate, feedback::Vector; skip_subspaces = [])
    time = [0.0]
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
    param[:t] = 80
    param[:dt] = 0.02
    param[:gamma] = 0.03
    param[:d] = 3
    n = 2
    param[:state] = ValueInfo(bosonstack(n, L, 1; d = 3), "|$(n)00...>")
    param[:alg] = "krylov BN-subspace"
    param[:k] = 6
    param[:H] = "BH d_w fixed |2> energy"
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "msr_n ends -> fb -> 0")
    return param
end

function f()
    

    # -----------
    Ls = [5]
    traj = 1
    # -----------
    plots = []
    param = initial_parameters(Ls[1])
    for L in Ls  
        param = initial_parameters(L)
        d = param[:d]
        indices = total_boson_number_subspace_indices(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        initial_id = find_subspace(param[:state].val, indices)
        state = subspace_split(param[:state].val, ranges, perm_mat)

        msrop = measurement_subspace(param[:effect].val.msr, ranges, perm_mat)
        feedback = feedback_measurement_subspace(param[:effect].val.fb, msrop, indices; digit_error = 10, id_relative_guess = -1)

        n_op = [subspace_split(singlesite_n(d, L, i), ranges, perm_mat) for i in 1:L]
        obs = [(state, id) -> expval(state[id], n_op[i][id]) for i in 1:L]

        dt = param[:dt]; t = param[:t]; k = param[:k]
        pa_k = PA_krylov_sub(state, k)
        pa_args, pa_out = pre_alloc_threads(zeros(length(obs), length(0:dt:t)), traj, pa_k)

        dw = param[:w].val .+ (rand(L) .- 0.5) * param[:dw].val
        dU = 2 * dw .- (2 * param[:w].val - param[:U].val)
        H = subspace_split(bosehubbard(d, L; w = dw, U = dU, J = 1), ranges, perm_mat)
        
        r_f(out, pa_k) = begin 
            effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            krylovevolve(state, initial_id, H, dt, t, k, pa_k, obs...; effect!, out)
        end

        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = traj_mean(pa_out)
        t_real = collect(0:dt:t) / (param[:J].val * 2 * π)
        map = heatmap(1:L, t_real, transpose(r), 
        xlabel = "sites", 
        ylabel = L"$\textrm{Time} \ (t) \ [\textrm{\mu s}]$",
        clim = (0, 2),
        #colorbar_ticks = ([0.0, 1.0, 2.0], ["0", "1", "2"]),
        #title = PlotAndSave.parameters_to_title(; L, param...)
        ) # blue #0472BE
        push!(plots, map)
    end

    # title = PlotAndSave.parameters_to_title(; param...)
    pl = plot(plots..., layout = (1, length(Ls)), size = (1500, 1500), margin = 20mm, 
    title = "Stack of excitations",
    titlecolor = :grey)
    savefig(pl, joinpath(@__DIR__, "heatmap.png"))
    savefig(pl, joinpath(@__DIR__, "heatmap.svg"))
    return pl
end

Plots.scalefontsizes(6)
@time pl = f()
Plots.scalefontsizes()
pl