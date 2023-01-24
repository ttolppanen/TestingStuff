using Distributed
@everywhere begin
    using Pkg; Pkg.activate(".")
    # Pkg.instantiate(); Pkg.precompile()
end
using QuantumStates
using QuantumOperators
using ITensors
using Gadfly
@everywhere using QuantumTimeEvolution

function make_title(; kwargs...)
    out = ""
    for (key, value) in kwargs
        out *= "$key = $value, "
    end
    return out[1:end-2]
end

function krylov(d, L, dt, t, k, traj, prob; color = :auto, kwargs...)
    function make_msrop(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    state0 = zeroone(d, L)
    H = bosehubbard(d, L; kwargs...)
    msrop = measurementoperators(nop(d), L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(d, L, s, Int(floor(L/2))))
    res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent; paral = :distributed)
    return layer(x = prob, y = res, Geom.line)
end

function mps(d, L, dt, t, traj, prob; addplot = true, color = :auto, kwargs...)
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
    addplot ? plot!(prob, res, label = "mps"; color) : plot(prob, res, label = "mps"; color)
end

function f()
    prob = 0.00:0.01:0.03
    traj = 100
    d = 2
    t = 30
    # @time krylov(d, 6, 0.02, t, 6, traj, prob; addplot = false, color = :red)
    @time l1 = krylov(d, 8, 0.02, t, 10, traj, prob; w=1)
    @time l2 = krylov(d, 8, 0.02, t, 10, traj, prob; w=10)
    # @time krylov(d, 10, 0.02, t, 6, traj, prob; addplot = true, color = :red)
    # @time krylov(d, 6, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time krylov(d, 8, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time krylov(d, 10, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time mps(2, 8, 0.6, 30, 1000, 0.01:0.1:0.91; addplot = true, maxdim = 5, cutoff = 1E-8)
    # title!(make_title(; d, t, traj) * ", std_msr")
    # figure_path = joinpath(dirname(@__FILE__), "same_gamma_small_dt_stantard_msr")
    # png(figure_path)
    plot(l1, l2)
end

f()
