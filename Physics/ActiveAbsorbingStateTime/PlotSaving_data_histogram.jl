using Distributed
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function return_times(result, dt; width = 0.01)
    out = Float64[]
    for traj in result
        did_not_return = true
        for i in 1:length(traj[1, :])
            time_step_value = traj[1, i]
            if abs(1 - time_step_value) <= width
                push!(out, i * dt)
                did_not_return = false
                break
            end
        end
        if did_not_return
            push!(out, -1)
        end
    end
    return out
end

function krylov(d, L, dt, t, traj, p; w, J, U, k)
    println("Generating")
    # state0 = onezero(d, L)
    state0 = 1im * normalize(ones(d^L))
    msrop = measurementoperators(nop(d), L)
    feedback = [isodd(i) ? n_bosons_projector(2, 0) : n_bosons_projector(2, 1) for i in 1:L]
    unitary_after_measurement!(feedback, msrop)
    bra_zeroone = zeroone(d, L)'
    calc_dst(state) = abs(bra_zeroone * state)
    H = bosehubbard(d, L; w, J, U)

    effect!(state) = measuresitesrandomly!(state, msrop, p)
    r_f() = krylovevolve(state0, H, dt, t, k, calc_dst; effect!)
    r = solvetrajectories(r_f, traj; paral = :distributed)
    return return_times(r, dt)
end

function f(traj)
    param = Dict(
        :dt => 0.02,
        :t => 120,
        :w => 0,
        :U => 1,
        :J => 1,
        :d => 2,
        :measurement => "std_n_feedback_01",
        :save_before_effect => false,
        :state => "total superposition",
        :alg => "krylov",
        :k => 6,
        :L => 12,
        :p => 0.06
    )
    data = krylov(param[:d], param[:L], param[:dt], param[:t], traj, param[:p]; w = param[:w], U = param[:U], J = param[:J], k = param[:k])
    println("Saving Data")
    path = joinpath(@__DIR__, "dst_histogram_STotalSuperposition_L12_$(param[:p])")
    makehistogram(path, data; traj, param...)
end

@time f(1000)
