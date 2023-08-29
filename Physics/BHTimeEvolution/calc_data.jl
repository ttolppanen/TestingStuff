using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave

function dis()
    return 2 * rand() - 1
end

function f()
    d = 3; L = 3;
    dt = 0.02; t = 10.0; k = 6;
    w = [0, 0, 0]; U = 0; J = 1;
    state0 = productstate(d, [0, 2, 0])
    H = bosehubbard(d, L, w; U, J)
    n1 = singlesite([0 0 0; 0 1 0; 0 0 0], L, 1)
    n2 = singlesite([0 0 0; 0 0 0; 0 0 1], L, 1)
    r = krylovevolve(state0, H, dt, t, k, s -> expval(s, n1), s -> expval(s, n2))
    
    x = 0:dt:t
    line1 = LineInfo(x, r[1, :], 1, "1_site_1_bosons")
    line2 = LineInfo(x, r[2, :], 1, "1_site_2_bosons")
    path = joinpath(@__DIR__, "first_site_boson_numbers_disorder")
    makeplot(path, line1, line2; d, L, dt, t, k, w = "disorder", U, J, state = "[020]")
end

f()