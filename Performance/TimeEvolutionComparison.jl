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
        if dif < 1E-8 break end
        ent1 .= ent2
    end
    return maxdim - 1
end

function f(alg)
    d = 3; L = 4
    dt = 0.02; t = 20
    repetitions = 10
    state0 = zeroone(d, L)
    mps0 = zeroonemps(d, L)
    H = bosehubbard(d, L)
    U = exp(-im * dt * Matrix(H))
    gates = bosehubbardgates(siteinds(mps0), dt)
    t = @elapsed begin
        for _ in 1:repetitions
            if alg == "e"
                exactevolve(state0, U, dt, t)
            elseif alg == "k"
                krylovevolve(state0, H, dt, t, 6)
            elseif alg == "m"
                mpsevolve(mps0, gates, dt, t)
            end
        end
    end
    println(alg * ": " * string(round(t/repetitions; digits=5)) * "s")
    # saveplot(pl, "middle_bsn_n")
end

# f("e")
# f("k")
# f("m")
# @show exactmpsdifference(3, 6;cutoff=1E-14, maxdim=10);
@show mpsconvergence(2, 12; maxdim = 32)