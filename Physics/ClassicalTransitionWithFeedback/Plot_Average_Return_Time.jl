using Distributed
using PlotAndSave
using Statistics
@everywhere using QuantumTimeEvolution

@everywhere begin 
function scramble_to_left!(state, i)
    if state[i - 1] == 0
        state[i] = 0
        state[i - 1] = 1
    end
end

function scramble_to_right!(state, i)
    if state[i + 1] == 0
        state[i] = 0
        state[i + 1] = 1
    end
end

function nearest_neighbor_scramble!(state, p)
    L = length(state)
    for i in 1:L
        if state[i] == 1 && rand() < p #if we have something to scramble, we scramble it with probability p.
            if i == 1
                scramble_to_right!(state, i)
            elseif i == L
                scramble_to_left!(state, i)
            else
                sites_to_scramble = [state[i - 1] == 0, state[i + 1] == 0]
                if all(sites_to_scramble)
                    if rand() < 0.5
                        scramble_to_right!(state, i)
                    else
                        scramble_to_left!(state, i)
                    end
                elseif sites_to_scramble[1]
                    scramble_to_left!(state, i)
                elseif sites_to_scramble[2]
                    scramble_to_right!(state, i)
                end
            end
        end
    end
end

function force_to_zeroone!(state, p)
    L = length(state)
    for i in 1:L
        if rand() < p # if we force.
            if isodd(i)
                state[i] = 0
            else
                state[i] = 1
            end
        end
    end
end

function classical_trajectory(initial_state, steps, scramble_prob, force_prob, observable)
    state = copy(initial_state)
    out = zeros(steps)
    for i in 1:steps
        nearest_neighbor_scramble!(state, scramble_prob)
        force_to_zeroone!(state, force_prob)
        out[i] = observable(state)
    end
    return out
end

function mean_return_time(result)
    out = []
    for traj in result
        did_not_return = true
        for i in eachindex(traj)
            if traj[i] == 1.0
                push!(out, i)
                did_not_return = false
                break
            end 
        end
        if did_not_return
            push!(out, length(traj))
        end
    end
    return mean(out)
end

end # @everywhere

function f()
    param = Dict(
        :steps => 10000,
        :scramble_prob => 0.9,
        :initial_state => "[101010...]",
        :traj => 1000
    )
    lines = []
    for L in 10:10:100
        initial_state = [isodd(i) ? 1 : 0 for i in 1:L]
        observable(state) = state == [isodd(i) ? 0 : 1 for i in 1:L]
        outcomes = []
        force_probs = 0.01:0.025:0.91
        for force_prob in force_probs
            r_f() = classical_trajectory(initial_state, param[:steps], param[:scramble_prob], force_prob, observable)
            r = solvetrajectories(r_f, param[:traj]; paral = :distributed)
            push!(outcomes, mean_return_time(r))
        end
        push!(lines, LineInfo(force_probs, outcomes, param[:traj], "L = $L"))
    end
    path = joinpath(@__DIR__, "large_systems")
    makeplot(path, lines...; xlabel = "force_prob", ylabel = "mean_steps", param...)
end

@time f()