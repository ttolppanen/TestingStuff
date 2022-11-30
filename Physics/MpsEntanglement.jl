using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots
using Random
include("Utility/PlotSaving.jl")

function f()
    d = 2; L = 8; half = Int(floor(L/2))
    dt = 0.1; t = 5.0; ts = 0:dt:t
    p=0
    state = zeroone(d, L)
    mps = onezeromps(d, L)
    H = bosehubbard(d, L)
    gates = bosehubbardgates(siteinds(mps), dt; k=1)
    effect! = measuresitesrandomly!(L, nop(d), p)
    meffect! = measuresitesrandomly!(siteinds(mps), nop(d), p)
    Random.seed!(2)
    r = krylovevolve(state, H, dt, t, 6; effect!)
    pl = plot(ts, entanglement(d, L, r, half), label="exact", title="msr_prob = $p")
    for maxdim in 4:4:16
        Random.seed!(2)
        r = mpsevolve(mps, gates, dt, t; maxdim, effect! = meffect!)
        plot!(pl, ts, entanglement(r, half), label="maxdim = $maxdim")
    end
    pl
    saveplot(pl, "entanglement_maxdim_p$p")
end

f()