using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using Plots

include("Utility/PlotSaving.jl")

function exactmpsdifference(d::Integer, L::Integer; kwargs...)
    half = Int(floor(L/2))
    dt = 0.01; t = 20
    state0 = zeroone(d, L)
    mps0 = zeroonemps(d, L)
    H = bosehubbard(d, L)
    U = exp(-im * dt * Matrix(H))
    gates = bosehubbardgates(siteinds(mps0), dt)
    r_e = krylovevolve(state0, H, dt, t, 6)
    r_m = mpsevolve(mps0, gates, dt, t; kwargs...)
    ent1 = [entanglement(d, L, s, half) for s in r_e]
    ent2 = [entanglement(s, half) for s in r_m]
    return maximum(abs.(ent1 .- ent2))
end
function mpsconvergence(d::Integer, L::Integer; maxdim=2)
    half = Int(floor(L/2))
    dt = 0.01; t = 20
    mps0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt)
    r_m = mpsevolve(mps0, gates, dt, t; maxdim)
    ent1 = [entanglement(s, half) for s in r_m]
    while true
        maxdim += 1
        @show maxdim
        r_m = mpsevolve(mps0, gates, dt, t; maxdim)
        ent2 = [entanglement(s, half) for s in r_m]
        dif = maximum(abs.(ent1 .- ent2))
        @show dif
        if dif < 1E-8 break end
        ent1 .= ent2
    end
    return maxdim - 1
end

function calctimes(L)
    d = 3
    dt = 0.02; t = 20
    # maxdim = Int((d^L)/2)
    maxdim = 10
    state0 = zeroone(d, L)
    mps0 = zeroonemps(d, L)
    H = bosehubbard(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt; k=2)
    t_k = @elapsed krylovevolve(state0, H, dt, t, 6)
    t_m = @elapsed mpsevolve(mps0, gates, 0.6, t; maxdim)
    return t_k, t_m
end

function mpsmaxdim(L, maxdim)
    d = 2
    dt = 0.6; t = 20
    mps0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt; k=2)
    t_m = @elapsed mpsevolve(mps0, gates, 0.6, t; maxdim)
    return t_m
end
function testingmaxdimdifference()
    m5_times = []
    m10_times = []
    m15_times = []
    Ls = 12:18
    for L in Ls
        t_m5 = mpsmaxdim(L, 5)
        t_m10 = mpsmaxdim(L, 10)
        t_m15 = mpsmaxdim(L, 15)
        push!(m5_times, t_m5)
        push!(m10_times, t_m10)
        push!(m15_times, t_m15)
    end
    pl = plot(Ls, m5_times, label = "maxdim = 5", yaxis=:log)
    plot!(pl, Ls, m10_times, label = "maxdim = 10")
    plot!(pl, Ls, m15_times, label = "maxdim = 15")
    # saveplot(pl, "mps maxdim")
    pl
end

function f()
    k_times = []
    m_times = []
    Ls = 4:12
    for L in Ls
        t_k, t_m = calctimes(L)
        push!(k_times, t_k)
        push!(m_times, t_m)
    end
    pl = plot(Ls, k_times, title = "d = 3, maxdim = 10", label = "krylov", yaxis=:log)
    plot!(pl, Ls, m_times, label = "mps")
    saveplot(pl, "krylov mps qutrit")
    pl
end

f()
# @show exactmpsdifference(2, 8; cutoff=1E-14, maxdim=16);
# @show mpsconvergence(2, 8; maxdim = 8)