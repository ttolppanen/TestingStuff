using Distributed
@everywhere using QuantumStates
@everywhere using QuantumOperators
@everywhere using ITensors
using PlotAndSave
@everywhere using QuantumTimeEvolution

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function krylov_with_disorder(d, L, dt, t, k, w, traj, prob; kwargs...)
    function mutual_information_half_centered(d, L)
        a = ceil(Int, L / 2) - floor(Int, L / 4) # ceil is for the odd-case
        b = L - (a - 1)
        return s -> entanglement(d, L, s, [a]) + entanglement(d, L, s, [b]) - entanglement(d, L, s, [a, b])
    end
    function mutual_information_first_to_half(d, L)
        half = Int(floor(L/2))
        return s -> entanglement(d, L, s, [1]) + entanglement(d, L, s, [half + 1]) - entanglement(d, L, s, [1, half + 1])
    end
    state0 = zeroone(d, L)
    msrop = measurementoperators(nop(d), L)
    calc_ent(state) = entanglement(d, L, state, Int(floor(L/2)))
    calc_mutual = mutual_information_first_to_half(d, L)
    make_H() = bosehubbard(d, L, w * (2 .* rand(L) .- 1); kwargs...)
    y = []
    for p in prob
        effect!(state) = measuresitesrandomly!(state, msrop, p)
        r_f() = krylovevolve(state0, make_H(), dt, t, k, calc_mutual; effect!)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        out = traj_mean(r)[1, end]
        push!(y, out)
    end
    return layer(x = prob ./ dt, y = y, Geom.line, color = ["dt = $dt, L = $L"])
end

function mps(d, L, dt, t, traj, prob; w, kwargs...)
    function mutual_information_first_to_half(d, L)
        half = Int(floor(L/2))
        return s -> entanglement(d, L, s, [1]) + entanglement(d, L, s, [half + 1]) - entanglement(d, L, s, [1, half + 1])
    end
    mps0 = zeroonemps(d, L)
    msrop = measurementoperators(nop(d), siteinds(mps0))
    calc_ent(state) = entanglement(state, Int(floor(L/2)))
    calc_mutual_f = mutual_information_first_to_half(d, L)
    calc_mutual(s) = calc_mutual_f(mps_to_array(s))
    make_gates() = bosehubbardgates(siteinds(mps0), dt, w * (2 .* rand(L) .- 1))

    mi_y = []
    ent_y = []
    for p in prob
        effect!(state) = measuresitesrandomly!(state, msrop, p; kwargs...)
        r_f() = mpsevolve(mps0, make_gates(), dt, t, calc_mutual, calc_ent; effect!, save_only_last = true, kwargs...)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        out = traj_mean(r)
        push!(mi_y, out[1, 1])
        push!(ent_y, out[2, 1])
    end
    return [LineInfo(prob ./ dt, mi_y, traj, "L = $L"), LineInfo(prob ./ dt, ent_y, traj, "L = $L")]
end

function f()
    dt = 0.3; t = 30
    prob = 0.03:0.06:0.63
    maxdim = 20
    cutoff = 1E-8
    disorder = 1
    d = 2
    traj = 10
    mi_layers = []
    ent_layers = []
    for L in 12:2:18
        @time l = mps(d, L, dt, t, traj, prob; w = disorder, maxdim, cutoff)
        push!(mi_layers, l[1])
        push!(ent_layers, l[2])
    end
    path = joinpath(@__DIR__, "MI_MPS")
    PlotAndSave.Plot(path, mi_layers...; xlabel = "gamma", ylabel = "MI(1, L/2)", d, dt, t, traj, maxdim, cutoff, disorder, save_before_msr = false)
    path = joinpath(@__DIR__, "ENT_MPS")
    PlotAndSave.Plot(path, ent_layers...; xlabel = "gamma", ylabel = "S", d, dt, t, traj, maxdim, cutoff, disorder, save_before_msr = false)
end

f()
