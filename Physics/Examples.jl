using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave
using Plots # <- julia package for plotting

function example()
    d = 3 # qubits
    L = 2 # number of sites
    dt = 0.4 * 50
    t = 1000 * 50
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = bosonstack(2, L, 1; d = 3) # |01> 
    w = 7500 / 250
    U = 1
    J = 1 / 250
    dw = w .+ (rand(L) .- 0.5) * 100 / 250
    dU = 2 * dw .- (2 * w - 1)
    H = bosehubbard(d, L; w=dw, U=dU, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    n = nall(d, L)
    calc_n_expval(state) = expval(state, n) # notice that calc_n_expval is a function with one argument, and n is "captured".
    n1 = singlesite_n(d, L, 1)
    calc_n1_expval(state) = expval(state, n1) # first site boson occupation number

    U = exp(-1.0im * Matrix(H) * dt)
    results = exactevolve(initial_state, U, dt, t, calc_n_expval, calc_n1_expval) # this returns a matrix, where the different columns refer to timesteps, and the rows to different observables.

    time_steps = 0:dt:t
    n_results = results[1, :]
    n1_results = results[2, :]
    
    pl = plot(time_steps, n_results, label = "Total boson number")
    plot!(pl, time_steps, n1_results, label = "First site boson number")

    # making the plot with PlotAndSave
    # plot_lines = [
    #     LineInfo(time_steps, n_results, 1, "Total boson number"), # here the number 1 refers to the number of trajectories.
    #     LineInfo(time_steps, n1_results, 1, "First site boson number") 
    # ]
    # save_data_path = joinpath(@__DIR__ , "plot_data")
    # makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "I dont know", d, L, dt, t, k, w, U, J, state = "|01>")
end

@time example()