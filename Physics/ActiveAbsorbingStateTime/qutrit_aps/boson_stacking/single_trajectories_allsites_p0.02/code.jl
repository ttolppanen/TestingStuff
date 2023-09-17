using Distributed
using JLD2
using Plots
using LinearAlgebra
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

@everywhere function measurement_and_feedback!(state::AbstractVector{<:Number}, msr_op, msr_prob::Real, fb_op)
    for i in 1:length(msr_op)
        if rand(Float64) < msr_prob
            measuresite!(state, msr_op, i)
            state .= fb_op[i] * state
            normalize!(state)
        end
    end
end

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function zeroone_orderparameter(L)
    d = 2 # not well defined for higher d
    order_parameter = singlesite_n(d, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        if isodd(i)
            order_parameter .+= singlesite_n(d, L, i)
        else
            order_parameter .+= singlesite([1 0; 0 1], L, 1) .- singlesite_n(d, L, i)
        end
    end
    order_parameter ./= L
    return state -> expval(state, order_parameter)
end

function n_squared_orderparameter(d, L)
    n2 = nop(d).^2
    order_parameter = singlesite(n2, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        order_parameter .+= singlesite(n2, L, i)
    end
    order_parameter ./= (2 * L)
    return state -> expval(state, order_parameter)
end

function anharmonicity_orderparameter(d, L)
    local_op = nop(d) * (nop(d) - I)
    out = singlesite(local_op, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        out .+= singlesite(local_op, L , i)
    end
    out ./= L
    return state -> expval(state, out)
end

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    state0 = zeroone(d, L)
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 0; 0 0 0; 0 0 1]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    n_list = [singlesite_n(d, L, i) for i in 1:L]
    calc_list = [s -> expval(s, n) for n in n_list]
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measurement_and_feedback!(state, msrop, p, fb_op)
    r = krylovevolve(state0, H, dt, t, k, calc_list...; effect!)
    return [LineInfo(0:dt:t, r[i, :], 1, "N_$i") for  i in 1:L]
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 10,
        :J => 1,
        :d => 3,
        :measurement => "std_n fb: 0 -> |0>, 1 -> |0>, 2 -> |2>",
        :save_before_effect => false,
        :state => "zeroone",
        :alg => "krylov",
        :L => 8,
        :k => 6,
        :p => 0.02
    )
    traj = 1
    order_line = krylov(param[:d], param[:L], param[:dt], param[:t], traj, param[:p]; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
    
    println("Saving Data")
    path = joinpath(@__DIR__, "boson_stacking" , "single_trajectories_allsites_p$(param[:p])")
    makeplot(path, order_line...; xlabel = "t", ylabel = "N_i", param...)
end


@time f()
