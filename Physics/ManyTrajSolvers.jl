using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function exact(d, L, dt, t, traj, msr_prob; addplot=false)
    state = zeroone(d, L)
    H = bosehubbard(d, L)
    U = exp(-im * dt * Matrix(H))
    effect! = measuresitesrandomly!(L, nop(d), msr_prob)
    r_f() = exactevolve(state, U, dt, t; effect!)
    r = solvetrajectories(r_f, traj)
    ent = trajmean(r, s -> entanglement(d, L, s, Int(floor(L/2))))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

function krylov(d, L, dt, t, k, traj, msr_prob; addplot=false)
    state = zeroone(d, L)
    # H = bosehubbard(d, L)
    effect! = measuresitesrandomly!(L, nop(d), msr_prob)
    r_f() = krylovevolve(state, bosehubbard(d, L, 2 * rand(L) .- 1), dt, t, k, s -> entanglement(d, L, s, Int(floor(L/2))); effect!)
    r = solvetrajectories(r_f, traj)
    ent = traj_mean(r)[1, :]
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

function make_msrop(d, state, L)
    msrop = measurementoperators(nop(d), siteinds(state))
    proj_op1 = n_bosons_projector(d, 1)
    proj_op2 = n_bosons_projector(d, 0)
    v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
    unitary_after_measurement!(v_p, msrop)
    return msrop
end

function mps(d, L, dt, t, traj, msr_prob; addplot=false, kwargs...)
    state = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(state), dt; k=2)
    msrop = make_msrop(d, state, L)
    effect!(state) = measuresitesrandomly!(state, msrop, msr_prob; kwargs...)
    # effect! = measuresitesrandomly!(siteinds(state), nop(d), msr_prob; maxdim)
    r_f() = mpsevolve(state, gates, dt, t; effect!, kwargs...)
    r = solvetrajectories(r_f, traj)
    ent = trajmean(r, s -> entanglement(s, Int(floor(L/2))))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

# exact(3, 3, 0.1, 5, 10, 0.1)
@time krylov(2, 6, 0.6, 30, 6, 100, 0.2; addplot=false)
# @time krylov(2, 6, 0.12, 30, 6, 2000, 0.11; addplot=true)
# @time mps(2, 6, 0.6, 30, 200, 0.05; addplot=true, maxdim=5, cutoff=1E-8)
