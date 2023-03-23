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

function mps(d, L, dt, t, traj, prob; w, J, U, kwargs...)
    function fluctuations(d, L)
        half = Int(floor(L/2))
        n = singlesite_n(d, L, 1)
        for i in 2:half
            n .+= singlesite_n(d, L, i)
        end
        n_squared = n * n
        f(s) = expval(s, n_squared) - expval(s, n)^2
        return s -> f(mps_to_array(s))
    end
    println("Generating")
    mps0 = zeroonemps(d, L)
    msrop = measurementoperators(nop(d), siteinds(mps0))
    calc_ent(state) = entanglement(state, Int(floor(L/2)))
    calc_fluc = fluctuations(d, L)
    gates = bosehubbardgates(siteinds(mps0), dt; w, U, J)

    ent_y = []
    fluc_y = []
    println("Going to Probabilities")
    for p in prob
        effect!(state) = measuresitesrandomly!(state, msrop, p; kwargs...)
        r_f() = mpsevolve(mps0, gates, dt, t, calc_ent, calc_fluc; effect!, save_before_effect = true, save_only_last = true, kwargs...)
        r = solvetrajectories(r_f, traj; paral = :distributed)
        out = traj_mean(r)
        push!(ent_y, out[1, 1])
        push!(fluc_y, out[2, 1])
    end
    return [LineInfo(prob, ent_y, traj, "L = $L"), LineInfo(prob, fluc_y, traj, "L = $L")]
end

function f(traj)
    param = Dict(
        :dt => 0.3,
        :t => 30,
        :cutoff => 1E-8,
        :maxdim => 90,
        :w => 0,
        :U => 5,
        :J => 1,
        :d => 3,
        :measurement => "std_n",
        :save_before_effect => true,
        :state => "zeroone",
        :trotter_order => 2,
        :alg => "MPS"
    )
    prob = 0.01:0.01:0.5
    fluc_layers = []
    ent_layers = []
    for L in [6, 10, 14]
        @time l = mps(param[:d], L, param[:dt], param[:t], traj, prob; w=param[:w], U=param[:U], J=param[:J], maxdim = param[:maxdim], cutoff = param[:cutoff])
        push!(ent_layers, l[1])
        push!(fluc_layers, l[2])
    end
    println("Saving Data")
    path = joinpath(@__DIR__, "qutrit_ent_90")
    addtrajectories(path, ent_layers...)
    # makeplot(path, ent_layers...; xlabel = "p", ylabel = "S", param...)
    path = joinpath(@__DIR__, "qutrit_fluc_90")
    addtrajectories(path, fluc_layers...)
    # makeplot(path, fluc_layers...; xlabel = "p", ylabel = "Fluc", param...)
end

@everywhere begin
    ITensors.set_warn_order(20)
end

@time f(200)
