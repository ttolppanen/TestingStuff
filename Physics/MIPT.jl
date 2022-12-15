using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots
include("Utility/PlotSaving.jl")

function krylov(d, L, dt, t, k, traj, prob; addplot = true)
    function make_msrop(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    state0 = zeroone(d, L)
    H = bosehubbard(d, L)
    msrop = measurementoperators(nop(d), L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(d, L, s, Int(floor(L/2))))
    res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent)
    addplot ? plot!(prob, res, label = "dt = $dt") : plot(prob, res, label = "dt = $dt")
end

function mps(d, L, dt, t, traj, prob; addplot = true, kwargs...)
    function make_msrop(d, state, L)
        msrop = measurementoperators(nop(d), siteinds(state))
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    mps0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt)
    msrop = measurementoperators(nop(d), siteinds(mps0))
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob; kwargs...)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(s, Int(floor(L/2))))
    res = mipt(mps0, gates, meffect!, dt, t, prob, traj, calc_ent; kwargs...)
    addplot ? plot!(prob, res, label = "L = $L") : plot(prob, res, label = "L = $L")
end


@time krylov(2, 6, 0.6, 30, 6, 1000, 0.1:0.1:0.4; addplot = false)
@time mps(2, 6, 0.6, 30, 100, 0.1:0.1:0.4; addplot = true, maxdim = 5, cutoff = 1E-8)