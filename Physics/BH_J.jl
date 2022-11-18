using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Plots

include("Utility/PlotSaving.jl")

function f()
    d = 3; L = 2
    dt = 0.01; t = 3
    J = 3
    target = 1
    state0 = singleone(d, L, target)
    result = exactevolve_bosehubbard(d, L, state0, dt, t; J)
    n = singlesite_n(d, L, target)
    result = expval(result, n)
    x_t = (0:dt:t) / (pi / J)
    pl = plot(x_t, result, title = "J = $(J)")
    #saveplot(pl, "middle_bsn_n")
    pl
end

f()