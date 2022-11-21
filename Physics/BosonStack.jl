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
    U = 50
    J = 2
    dt = 0.1 * U; t = 15 * U
    J_e = (-1)^(N - 1) * N / factorial(N - 1) * (J/U)^(N - 1) * J
    target = 2
    state0 = bosonstack(N, L, target)
    result = exactevolve_bosehubbard(d, L, state0, dt, t; U, J)
    res = []
    for i in 1:4
        push!(res, n_i_expval(result, d, L, i))
    end
    times = (0:dt:t) / (pi / J_e)
    pl = plot(times, res[1], label = "site = 1", title = "U = $(U), J = $(J)")
    plot!(pl, times, res[2], label = "site = 2")
    plot!(pl, times, res[3], label = "site = 3")
    plot!(pl, times, res[4], label = "site = 4")
    plot!(pl, times, x -> 3 * (cos(x * pi)^2))
    #saveplot(pl, "boson_stack_$(target)")
    pl
end

f()