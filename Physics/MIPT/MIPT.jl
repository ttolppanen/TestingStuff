using Distributed
@everywhere begin
    using Pkg; Pkg.activate(".")
    # Pkg.instantiate(); Pkg.precompile()
end
using QuantumStates
using QuantumOperators
using ITensors
using Plots
@everywhere using QuantumTimeEvolution

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
    msrop = make_msrop(d, L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(d, L, s, Int(floor(L/2))))
    res = mipt(state0, H, k, meffect!, dt, t, prob * dt, traj, calc_ent; paral = :distributed)
    color = dt == 0.3 ? :red : :blue 
    addplot ? plot!(prob, res, label = "L = $L, dt = $dt"; color) : plot(prob, res, label = "L = $L, dt = $dt"; color)
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
    res = mipt(mps0, gates, meffect!, dt, t, prob, traj, calc_ent; paral = :distributed, kwargs...)
    addplot ? plot!(prob, res, label = "mps") : plot(prob, res, label = "mps")
end

function f()
    prob = 1.0:0.01:1.61
    traj = 10000
    @time krylov(2, 6, 0.3, 30, 6, traj, prob; addplot = false)
    @time krylov(2, 8, 0.3, 30, 6, traj, prob; addplot = true)
    @time krylov(2, 10, 0.3, 30, 6, traj, prob; addplot = true)
    @time krylov(2, 6, 0.6, 30, 6, traj, prob; addplot = true)
    @time krylov(2, 8, 0.6, 30, 6, traj, prob; addplot = true)
    @time krylov(2, 10, 0.6, 30, 6, traj, prob; addplot = true)
    # @time mps(2, 8, 0.6, 30, 1000, 0.01:0.1:0.91; addplot = true, maxdim = 5, cutoff = 1E-8)
    figure_path = joinpath(dirname(@__FILE__), "same_gamma_collapse_zoom")
    png(figure_path)
end

f()
