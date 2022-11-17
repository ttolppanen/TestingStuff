using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Plots

include("Utility/PlotSaving.jl")

function n_i_expval(result, d, L, i)
    n = singlesite_n(d, L, i)
    return expval(result, n)
end

function f()
    d = 4; L = 4; N = 3
    dt = 0.01; t = 30
    U = 10
    J = 2
    target = 1
    state0 = bosonstack(N, L, target)
    result = exactevolve_bosehubbard(d, L, state0, dt, t; U, J)
    res = []
    for i in 1:4
        push!(res, n_i_expval(result, d, L, i))
    end
    pl = plot(0:dt:t, res[1], label = "site = 1", title = "U = $(U), J = $(J)")
    plot!(pl, 0:dt:t, res[2], label = "site = 2")
    plot!(pl, 0:dt:t, res[3], label = "site = 3")
    plot!(pl, 0:dt:t, res[4], label = "site = 4")
    saveplot(pl, "boson_stack_$(target)")
    pl
end

f()