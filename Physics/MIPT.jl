using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Plots
include("Utility/PlotSaving.jl")

function mipt(d, L, dt, t, traj, prob)
    state0 = zeroone(d, L)
    op_to_msr = nop(d)
    msrop = measurementoperators(op_to_msr, L)
    res = []
    for p in prob
        meffect!(state) = measuresitesrandomly!(state, msrop, p)
        r_exact() = krylovevolve_bosehubbard(d, L, state0, dt, t, 6; effect! = meffect!, savelast = true)
        r_traj = solvetrajectories(r_exact, traj)
        r_mean = trajmean(r_traj, s->entanglement(d, L, s, Int(floor(L/2))))
        push!(res, r_mean[1])
    end
    return prob, res
end

function f()
    res(L) = mipt(2, L, 0.1, 20.0, 50, 0.01:0.01:0.3)

    pl = plot(res(4))
    plot!(pl, res(6))
    plot!(pl, res(8))
    #saveplot(pl, "MIPT_qubits")
    pl
end

@time f()

