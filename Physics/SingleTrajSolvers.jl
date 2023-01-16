using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots

function exact(d, L, dt, t; addplot=false)
    state = zeroone(d, L)
    H = bosehubbard(d, L)
    U = exp(-im * dt * Matrix(H))
    r = exactevolve(state, U, dt, t)
    ent = entanglement(d, L, r, Int(floor(L/2)))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

function krylov(d, L, dt, t, k; addplot=false)
    state = zeroone(d, L)
    H = bosehubbard(d, L)
    r = krylovevolve(state, H, dt, t, k)
    ent = entanglement(d, L, r, Int(floor(L/2)))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

function mps(d, L, dt, t; addplot=false)
    state = zeroonemps(d, L)
    linkinds(state)
    gates = bosehubbardgates(siteinds(state), dt; k=2)
    r = mpsevolve(state, gates, dt, t; cutoff = 1E-8, maxdim = 40)
    ent = entanglement(r, Int(floor(L/2)))
    ts = 0:dt:t
    return addplot ? plot!(ts, ent) : plot(ts, ent)
end

# exact(3, 3, 0.1, 2)
# krylov(2, 8, 0.02, 20, 6; addplot=false)
@time mps(3, 15, 0.025, 10; addplot=false)
