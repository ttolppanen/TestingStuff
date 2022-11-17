using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Plots

include("Utility/PlotSaving.jl")

function f()
    d = 3; L = 5
    dt = 0.01; t = L
    J = 1
    target = 3
    state0 = singleone(d, L, target)
    result = exactevolve_bosehubbard(d, L, state0, dt, t; J)
    n = singlesite_n(d, L, target)
    result = expval(result, n)
    pl = plot(0:dt:t, result, title = "J = $(J)")
    saveplot(pl, "middle_bsn_n")
    pl
end

f()