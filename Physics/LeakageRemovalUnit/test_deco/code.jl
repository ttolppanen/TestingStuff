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
    ops = [sqrt(10) * singlesite_n(d, L, 1)]
    for op in ops
        H .-= 1.0im / 2 .* op' * op
    end
    sigma_x = sparse(complex([0 1.0; 1 0]))
    #sigma_x = singlesite(sigma_x, L, 1)
    obs = (state) -> expval(state, sigma_x)

    r_f() = krylovevolve(state, H, ops, dt, t, k, obs)
    r = solvetrajectories(r_f, 1000)

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

    push!(lines, LineInfo(x_t, r, 100, "data"))
    
    x_max = findlocalmaxima(r)
    push!(lines, LineInfo(x_t[x_max], r[x_max], 100, "maxima"))
    push!(lines, LineInfo(x_t, exp.(-x_t * 10), 1, "theoretical"))
    
    makeplot(joinpath(@__DIR__), lines...; copy_code = false)
end

f()