using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots
include("Utility/PlotSaving.jl")

function mipt(d, L, dt, t, traj, prob)
    cutoff=1E-6; maxdim=10
    state0 = zeroonemps(d, L)
    op_to_msr = nop(d)
    msrop = measurementoperators(op_to_msr, siteinds(state0))
    res = []
    for p in prob
        meffect!(state) = measuresitesrandomly!(state, msrop, p; cutoff, maxdim)
        r_f() = mpsevolve_bosehubbard(state0, dt, t; effect! = meffect!, savelast = true, cutoff, maxdim)
        r_traj = solvetrajectories(r_f, traj)
        r_mean = trajmean(r_traj, s->entanglement(s, Int(floor(L/2))))
        push!(res, r_mean[1])
    end
    return prob, res
end

function f()
    res(L) = mipt(2, L, 0.1, 20.0, 20, 0.01:0.1:0.3)

    pl = plot(res(8))
    #plot!(pl, res(6))
    #plot!(pl, res(8))
    #saveplot(pl, "MIPT_qubits")
    pl
end

@time f()

