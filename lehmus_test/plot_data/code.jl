using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave
using Plots # <- julia package for plotting

function example()
    d = 2 # qubits
    L = 2 # number of sites
    dt = 0.01
    t = 5.0
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = zeroone(d, L) # |01> 
    w = 0.0
    U = 2.0
    J = 3.0
    H = bosehubbard(d, L; w=w, U=U, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    msr_prob = 0.01
    effect! = random_measurement_function(L, nop(d), msr_prob)

    n1 = singlesite_n(d, L, 1)
    calc_n1_expval(state) = expval(state, n1) # first site boson occupation number

    results = krylovevolve(initial_state, H, dt, t, k, calc_n1_expval; effect!) # this returns a matrix, where the different columns refer to timesteps, and the rows to different observables.

    time_steps = 0:dt:t
    n1_results = results[1, :]

    plot_lines = [
        LineInfo(time_steps, n1_results, 1, "First site boson number") 
    ]
    save_data_path = joinpath(@__DIR__ , "plot_data")
    makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "I dont know", d, L, dt, t, k, w, U, J, state = "|01>")
end

@time example()