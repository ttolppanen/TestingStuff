using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots
include("Utility/PlotSaving.jl")

function mipt(d, L, dt, t, traj, prob)
    maxdim = 4
    state0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(state0), dt)
    op_to_msr = nop(d)
    msrop = measurementoperators(op_to_msr, siteinds(state0))
    res = []
    for p in prob
        meffect!(state) = measuresitesrandomly!(state, msrop, p; cutoff, maxdim)
        r_f() = mpsevolve(state0, gates, dt, t; effect! = meffect!, savelast = true, maxdim)
        r_traj = solvetrajectories(r_f, traj)
        r_mean = trajmean(r_traj, s->entanglement(s, Int(floor(L/2))))
        push!(res, r_mean[1])
    end
    return prob, res
end

function f()
    res(L) = mipt(2, L, 0.6, 20.0, 10, 0.2:0.2:0.8)

    pl = plot(res(8))
    plot!(pl, res(10))
    plot!(pl, res(12))
    # plot!(pl, res(10))
    # saveplot(pl, "MIPT_qubits")
    pl
end

# 30min for 200 traj, 8 10 12 qubits, 0:0.6:20.0t, 0.2:0.1:0.6
@profview @time f()

