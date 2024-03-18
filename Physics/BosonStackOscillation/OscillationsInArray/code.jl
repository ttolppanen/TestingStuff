using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave
using Plots # <- julia package for plotting
using NumericalIntegration

function example()
    d = 3 # qutrits
    L = 2 # number of sites
    dt = 0.01
    t = 20.0
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = productstate(d, [isodd(i) ? 2 : 0 for i in 1:L]) # |01> 
    w = 0.0
    U = 10.0
    J = 1.0
    H = bosehubbard(d, L; w=w, U=U, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    n_i = [singlesite_n(d, L, i) for i in 1:L]
    n = n_i[1] * (n_i[1] - I)
    for i in 2:L
        n .+= n_i[i] * (n_i[i] - I)
    end
    obs(state) = expval(state, n) / L

    n1 = singlesite_n(d, L, 1)
    n1_1bos = n1 * (n1 - I * 2)
    calc_n_expval(state) = expval(state, n1_1bos) # notice that calc_n_expval is a function with one argument, and n is "captured".

    results = krylovevolve(initial_state, H, dt, t, k, obs) # this returns a matrix, where the different columns refer to timesteps, and the rows to different observables.

    time_steps = 0:dt:t
    r = results[1, :]
    
    plot!(time_steps, r, label = "U = $U")

    # making the plot with PlotAndSave
    # plot_lines = [
    #     LineInfo(time_steps, n_results, 1, "Total boson number"), # here the number 1 refers to the number of trajectories.
    #     LineInfo(time_steps, n1_results, 1, "First site boson number") 
    # ]
    # save_data_path = joinpath(@__DIR__ , "plot_data")
    # makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "I dont know", d, L, dt, t, k, w, U, J, state = "|01>")
end

function add_vline_gamma(p)
    gamma = p / 0.02
    vline!([1 / gamma])
end

@time example()