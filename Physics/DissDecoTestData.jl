using DifferentialEquations
using QuantumStates
using QuantumOperators
using JLD2
using Plots
using LinearAlgebra

function dens_expval(rho, op)
    return real(tr(rho * op))
end

function f()
    d = 3; L = 3
    dt = 0.1; t = 100.0
    w = 2; U = 13

    s0 = bosonstack(2, 3, 1)
    rho0 = s0 * s0'
    H = bosehubbard(d, L; w, U)

    gamma = 0.004
    kappa = 0.02
    dd_op = vcat([sqrt(gamma) * singlesite_a(d, L, i) for i in 1:L], [sqrt(kappa) * singlesite_n(d, L, i) for i in 1:L])

    param = (H = H, dd_op = dd_op)

    prob = ODEProblem(diff_eq, rho0, (0.0, t), param; saveat = dt)
    sol = solve(prob,Tsit5())

    n1 = singlesite_n(d, L , 1)
    n1_res = [dens_expval(sol[:, :, i], n1) for i in 1:1:Int(t/dt)]
    n2 = singlesite_n(d, L , 2)
    n2_res = [dens_expval(sol[:, :, i], n2) for i in 1:1:Int(t/dt)]
    n3 = singlesite_n(d, L , 3)
    n3_res = [dens_expval(sol[:, :, i], n3) for i in 1:1:Int(t/dt)]
    n_tot = nall(d, L)
    n_tot_res = [dens_expval(sol[:, :, i], n_tot) for i in 1:1:Int(t/dt)]
    x = range(0.0, t, length = Int(t / dt))
    pl = plot(x, n1_res, label = "n1")
    plot!(pl, x, n2_res, label = "n2")
    plot!(pl, x, n3_res, label = "n3")
    plot!(pl, x, n_tot_res, label = "n tot")

    data = (t = x, n1 = n1_res, n2 = n2_res, n3 = n3_res, ntot = n_tot_res)
    jldsave("diss_deco_test_data"; data = data)
end

function diff_eq(drho, rho, param, t)
    drho .= -1im .* (param.H * rho - rho * param.H)
    for op in param.dd_op
        drho .+= 0.5 * (2 * op * rho * op' - op' * op * rho - rho * op' * op)
    end
end

sol = f()