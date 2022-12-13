using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Plots
include("Utility/PlotSaving.jl")

function make_msrop(d, L)
    msrop = measurementoperators(nop(d), L)
    proj_op1 = n_bosons_projector(d, 1)
    proj_op2 = n_bosons_projector(d, 0)
    v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
    # unitary_after_measurement!(v_p, msrop)
    return msrop
end

function f(d, L, dt, t, k, traj; addplot = true)
    gamma = 1:1:10
    prob = dt * gamma
    state0 = zeroone(d, L)
    H = bosehubbard(d, L)
    msrop = make_msrop(d, L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(d, L, s, Int(floor(L/2))))
    res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent)
    addplot ? plot!(gamma, res, label = "dt = $dt") : plot(gamma, res, label = "dt = $dt")
end

@time f(2, 6, 0.02, 30, 6, 5000; addplot = false)
@time f(2, 6, 0.04, 30, 6, 5000; addplot = true)
@time pl = f(2, 6, 0.08, 30, 6, 5000; addplot = true)
saveplot(pl, "mipt_timestep_difference")