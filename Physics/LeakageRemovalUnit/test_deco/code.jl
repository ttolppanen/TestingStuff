using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using LinearAlgebra
using PlotAndSave
using DSP
using SparseArrays
using Images

function f()
    d = 2; L = 1;
    dt = 0.02; t = 10;
    k = 6
    state = complex(sparse([1.0, 1]))
    normalize!(state)
    H = bosehubbard(d, L)
    deco_rate = 0.8
    ops = [sqrt(2 * deco_rate) * singlesite_n(d, L, 1)]
    for op in ops
        H .-= 1.0im / 2 .* op' * op
    end
    sigma_x = sparse(complex([0 1.0; 1 0]))
    #sigma_x = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, sigma_x)

    r_f() = krylovevolve(state, H, ops, dt, t, k, obs)
    r = solvetrajectories(r_f, 4000; paral = :threads)

    function traj_mean(result)
        out = zeros(size(result[1]))
        for traj in result
            out .+= traj
        end
        return out ./ length(result)
    end
    
    # r = map(x -> abs.(hilbert(x[1,:])), r)
    r = traj_mean(r)
    r = r[1, :]
    x_t = 0:dt:t
    lines = []

    pl = plot(x_t, r, label = "data")
    #push!(lines, LineInfo(x_t, r, 100, "data"))
    
    # x_max = findlocalmaxima(r)
    # push!(lines, LineInfo(x_t[x_max], r[x_max], 100, "maxima"))
    y_t(x) = cos(x) * exp(-x * deco_rate)
    #push!(lines, LineInfo(x_t, y_t.(x_t), 1, "theoretical"))
    plot!(pl, x_t, y_t, label = "analytical", linestyle = :dash)
    #makeplot(joinpath(@__DIR__), lines...; copy_code = false)
end

f()