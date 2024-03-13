using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave
using Plots
using NumericalIntegration

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function example()
    d = 3
    L = 6
    dt = 0.01
    t = 30
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = [0 for i in 1:L]
    initial_state[1] = 2
    initial_state = productstate(d, initial_state) 
    w = 0.0
    U = 15.0
    J = 1.0
    H = bosehubbard(d, L; w=w, U=U, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    n_L = singlesite_n(d, L, 1)
    calc_last_site_expval(state) = expval(state, n_L) # notice that calc_n_expval is a function with one argument, and n is "captured".
    msrop = measurementoperators(nop(d), L)
    effect!(state) = begin 
        measuresite!(state, msrop, 3)
        return state
    end

    r_f() = krylovevolve(initial_state, H, dt, t, k, calc_last_site_expval; effect!) # this returns a matrix, where the different columns refer to timesteps, and the rows to different observables.
    results = solvetrajectories(r_f, 100)
    results = mean(results)

    time_steps = 0:dt:t
    r = results[1, :]
    
    plot!(time_steps, r, label = "U = $U", ylabel = "N_$L")

    # making the plot with PlotAndSave
    # plot_lines = [
    #     LineInfo(time_steps, n_results, 1, "Total boson number"), # here the number 1 refers to the number of trajectories.
    #     LineInfo(time_steps, n1_results, 1, "First site boson number") 
    # ]
    # save_data_path = joinpath(@__DIR__ , "plot_data")
    # makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "I dont know", d, L, dt, t, k, w, U, J, state = "|01>")
end

@time example()