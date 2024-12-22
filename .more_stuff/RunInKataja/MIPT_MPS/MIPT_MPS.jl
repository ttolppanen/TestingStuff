using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots

function make_msrop(d, state, L)
    msrop = measurementoperators(nop(d), siteinds(state))
    proj_op1 = n_bosons_projector(d, 1)
    proj_op2 = n_bosons_projector(d, 0)
    v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
    unitary_after_measurement!(v_p, msrop)
    return msrop
end

function f(L::Integer; addplot = true)
    d = 2;
    dt = 0.6; t = 20
    prob = 0.1:0.1:0.9; traj = 1000
    cutoff = 1E-8; maxdim = 5
    mps0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt; U=1)
    msrop = make_msrop(d, mps0, L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob; cutoff, maxdim)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(s, Int(floor(L/2))))
    res = mipt(mps0, gates, meffect!, dt, t, prob, traj, calc_ent; cutoff, maxdim)
    addplot ? plot!(prob, res, label = "L = $L") : plot(prob, res, label = "L = $L")
end

@time f(6; addplot = false)
# @time f(8)
@time display(f(10))
savefig("mipt")