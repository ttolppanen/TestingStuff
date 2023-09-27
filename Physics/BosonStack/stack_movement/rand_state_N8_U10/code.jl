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

function anharmonicity_orderparameter(d, L)
    local_op = nop(d) * (nop(d) - I)
    out = singlesite(local_op, L , 1) #this is probably VERY slow for high d^L!
    for i in 2:L
        out .+= singlesite(local_op, L , i)
    end
    out ./= L
    return state -> expval(state, out)
end

@everywhere function rand_boson_state(d, L, N)
    state_fock = zeros(Int, L)
    while sum(state_fock) < N
        site = rand(1:L)
        if state_fock[site] < (d - 1)
            state_fock[site] += 1
        end
    end
    return productstate(d, state_fock)
end

function krylov(d, L, N, dt, t, traj; w, J, U, k)
    println("Generating")
    calc_order = anharmonicity_orderparameter(d, L)
    H = bosehubbard(d, L; w, J, U)

    r_f() = krylovevolve(rand_boson_state(d, L, N), H, dt, t, k, calc_order)
    println("Solving trajectories")
    r = solvetrajectories(r_f, traj; paral = :distributed)
    r = traj_mean(r)
    dst_y = r[1, :]
    return LineInfo(0:dt:t, dst_y, traj, "")
end

function f()
    param = Dict(
        :dt => 0.02,
        :t => 30,
        :w => 0,
        :U => 10,
        :J => 1,
        :d => 3,
        :save_before_effect => false,
        :state => "rand_N_boson_state",
        :alg => "krylov",
        :L => 8,
        :N => 8,
        :k => 6
    )
    traj = 100
    dst_lines = []
    order_line = krylov(param[:d], param[:L], param[:N], param[:dt], param[:t], traj; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
    push!(dst_lines, order_line)
    
    println("Saving Data")
    path = joinpath(@__DIR__, "rand_state_N$(param[:N])_U$(param[:U])")
    # addtrajectories(path, dst_lines...)
    makeplot(path, dst_lines...; xlabel = "t", ylabel = "sum(N(N-1)) / L", param...)
end

@time f()