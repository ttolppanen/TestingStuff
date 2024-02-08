using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end


function example()
    d = 2 # qubits
    L = 2 # number of sites
    dt = 0.01
    t = 2 * pi
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = onezero(d, L) # |01> 
    w = 1.0
    U = 0.0
    J = 1.0
    H = bosehubbard(d, L; w=w, U=U, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    n = nall(d, L)
    calc_n_expval(state) = expval(state, n) # notice that calc_n_expval is a function with one argument, and n is "captured".
    n1 = singlesite_n(d, L, 1)
    calc_n1_expval(state) = expval(state, n1) # first site boson occupation number
    calc_ent(state) = entanglement(d, L, state, 1)
    msrop = measurementoperators(nop(d), L)
    effect!(state) = random_measurement!(state, msrop, 0.005)

    r_f() = krylovevolve(initial_state, H, dt, t, k, calc_n1_expval, calc_ent; effect!) # this returns a matrix, where the different columns refer to timesteps, and the rows to different observables.

    results = r_f()
    time_steps = 0:dt:t
    # n_results = results[1, :]
    n1_results = results[1, :]
    ent_results = results[2, :]

    res_average = solvetrajectories(r_f, 100)
    res_average = traj_mean(res_average)
    n1_avg_results = res_average[1, :]

    plot_lines = [
        # LineInfo(time_steps, n_results, 1, "Total boson number"), # here the number 1 refers to the number of trajectories.
        LineInfo(time_steps, n1_results, 1, "First site excitation number")
        # LineInfo(time_steps, n1_avg_results, 1, "Average")
        # LineInfo(time_steps, ent_results, 1, "Entanglement")
    ]
    save_data_path = joinpath(@__DIR__ , "plot_data")
    makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "", d, L, dt, t, k, w, U, J, state = "|10>")
end

@time example()