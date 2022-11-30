using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using ITensors
using BenchmarkTools
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
    d = 2
    dt = 0.02; t = 20
    # maxdim = Int((d^L)/2)
    maxdim = 5
    cutoff = 1E-8
    state0 = zeroone(d, L)
    mps0 = zeroonemps(d, L)
    H = bosehubbard(d, L)
    U = exp(-im * dt * Matrix(H))
    gates = bosehubbardgates(siteinds(mps0), dt; k=1)
    t_e = @belapsed exactevolve($state0, $U, $dt, $t)
    t_k = @belapsed krylovevolve($state0, $H, $dt, $t, 6)
    t_m = @belapsed mpsevolve($mps0, $gates, $dt, $t; $cutoff, $maxdim)
    return t_e, t_k, t_m
end

function f()
    e_times = []
    k_times = []
    m_times = []
    Ls = 10:12
    for L in Ls
        t_e, t_k, t_m = calctimes(L)
        push!(e_times, t_e)
        push!(k_times, t_k)
        push!(m_times, t_m)
    end
    pl = plot(Ls, e_times, title = "maxdim = 5, cutoff = 1E-8", label = "exact")
    plot!(pl, Ls, k_times, label = "krylov")
    plot!(pl, Ls, m_times, label = "mps")
    # saveplot(pl, "time comparison")
end

f()
# @show exactmpsdifference(2, 8; cutoff=1E-14, maxdim=16);
# @show mpsconvergence(2, 8; maxdim = 8)