using Distributed
@everywhere begin
    using Pkg; Pkg.activate(".")
    # Pkg.instantiate(); Pkg.precompile()
end
using QuantumStates
@everywhere using QuantumOperators
using ITensors
using Gadfly
@everywhere using QuantumTimeEvolution

function make_title(; kwargs...)
    out = ""
    for (key, value) in kwargs
        out *= "$key = $value, "
    end
    return out[1:end-2]
end

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function krylov(d, L, dt, t, k, traj, prob; color = :auto, kwargs...)
    function make_msrop(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    state0 = zeroone(d, L)
    H = bosehubbard(d, L; kwargs...)
    msrop = measurementoperators(nop(d), L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(state) = entanglement(d, L, state, Int(floor(L/2)))
    res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent; paral = :distributed)
    return layer(x = prob ./ dt, y = res, Geom.line, color = ["dt = $dt, L = $L"])
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

function mps(d, L, dt, t, traj, prob; addplot = true, color = :auto, kwargs...)
    function make_msrop(d, state, L)
        msrop = measurementoperators(nop(d), siteinds(state))
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    mps0 = zeroonemps(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt)
    msrop = measurementoperators(nop(d), siteinds(mps0))
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob; kwargs...)
    calc_ent(state) = entanglement(state, Int(floor(L/2)))
    res = mipt(mps0, gates, meffect!, dt, t, prob, traj, calc_ent; paral = :distributed, save_before_effect = true, kwargs...)
    return layer(x = prob, y = res, Geom.line, color = ["L = $L"])
end

function f()
    dt = 0.16
    prob002 = 0.008:0.008:0.248
    traj = 2000
    d = 2
    t = 30
    disorder = 1
    # @time krylov(d, 6, 0.02, t, 6, traj, prob; addplot = false, color = :red)
    @time l1 = krylov_with_disorder(d, 6, dt, t, 10, disorder, traj, prob002)
    @time l2 = krylov_with_disorder(d, 8, dt, t, 10, disorder, traj, prob002)
    @time l3 = krylov_with_disorder(d, 10, dt, t, 10, disorder, traj, prob002)
    # @time l2 = krylov(d, 6, 0.32, t, 12, traj, prob032; w=0)
    # @time l11 = krylov(d, 8, 0.02, t, 6, traj, prob002; w=0)
    # @time l22 = krylov(d, 8, 0.32, t, 12, traj, prob032; w=0)
    # @time krylov(d, 10, 0.02, t, 6, traj, prob; addplot = true, color = :red)
    # @time krylov(d, 6, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time krylov(d, 8, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time krylov(d, 10, 0.2, t, 6, traj, prob; addplot = true, color = :blue)
    # @time mps(2, 8, 0.6, 30, 1000, 0.01:0.1:0.91; addplot = true, maxdim = 5, cutoff = 1E-8)
    # title!(make_title(; d, t, traj) * ", std_msr")
    # figure_path = joinpath(dirname(@__FILE__), "same_gamma_small_dt_stantard_msr")
    # png(figure_path)
    title = (make_title(; d, t, traj, k = 6) * ", std_msr, measurement before time-step")
    p = plot(l1, l2, l3, Guide.title(title), Guide.xlabel("Gamma"))
    p |> SVG("d_dt016_mifh.svg")
end

function calc_with_mps()
    dt = 0.3; t = 20
    prob = 0.15:0.075:0.9
    maxdim = 10
    cutoff = 1E-8
    d = 2
    traj = 500
    layers = []
    for L in 8:2:12
        @time l = mps(d, L, dt, t, traj, prob; maxdim, cutoff)
        push!(layers, l)
    end
    title = (make_title(; d, dt, t, traj, maxdim, cutoff) * ", \n std_msr, saving data before msr")
    p = plot(layers..., Guide.title(title), Guide.xlabel("p"))
    p |> SVG("test.svg", 6inch, 4inch)
end

f()
