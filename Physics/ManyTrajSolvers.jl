using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots

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
    H = bosehubbard(d, L)
    effect! = measuresitesrandomly!(L, nop(d), msr_prob)
    r_f() = krylovevolve(state, H, dt, t, k; effect!)
    r = solvetrajectories(r_f, traj)
    ent = trajmean(r, s -> entanglement(d, L, s, Int(floor(L/2))))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

function mps(d, L, dt, t, traj, msr_prob; addplot=false)
    maxdim = 4
    state = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(state), dt; k=2)
    effect! = measuresitesrandomly!(siteinds(state), nop(d), msr_prob; maxdim)
    r_f() = mpsevolve(state, gates, dt, t; effect!, maxdim)
    r = solvetrajectories(r_f, traj)
    ent = trajmean(r, s -> entanglement(s, Int(floor(L/2))))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

# exact(3, 3, 0.1, 5, 10, 0.1)
# @time krylov(2, 8, 0.02, 20, 6, 100, 0.02; addplot=false)
@time mps(2, 16, 0.6, 20, 100, 0.02 * 0.6 / 0.02; addplot=true)
