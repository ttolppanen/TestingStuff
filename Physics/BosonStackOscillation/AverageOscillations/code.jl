using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using PlotAndSave
using Plots # <- julia package for plotting
using NumericalIntegration

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function example()
    traj = 5000
    d = 3 # qutrits
    L = 2 # number of sites
    dt = 0.01
    t = 1.0
    k = 6 # the krylov subspace dimenion, this affects the accuracy of the numerics.
    initial_state = productstate(d, [isodd(i) ? 2 : 0 for i in 1:L]) # |01> 
    w = 0.0
    U = 7
    J = 1.0
    H = bosehubbard(d, L; w=w, U=U, J=J) # in julia the arguments after ; are keyword arguments. Here since the argument name is the same as the keyword, you dont need to write w=w, you can just put w. 

    one_one = allone(d, L)
    one_one = one_one * one_one'
    calc_one_one_expval(state) = expval(state, one_one) # notice that calc_n_expval is a function with one argument, and n is "captured".
    calc_ent(state) = entanglement(d, L, state, 1)

    p = 0.05
    op = nop(d)
    msr_op = measurementoperators(op, L)
    effect!(state) = random_measurement!(state, msr_op, p)

    r_f() = krylovevolve(initial_state, H, dt, t, k, calc_one_one_expval, calc_ent; effect!)
    results = solvetrajectories(r_f, traj)
    results = traj_mean(results)

    time_steps = 0:dt:t
    r = results[1, :]
    r_ent = results[2, :]
    
    int = float(integrate(time_steps, r))
    int = int / t
    int = round(int, digits = 3)
    
    #plot!(time_steps, r, label = "|11> density, U = $U, I/t = $int")
    plot!(time_steps, r_ent, label = "Ent, U = $U")
    #check_rabi(U, J, dt)
    # making the plot with PlotAndSave
    # plot_lines = [
    #     LineInfo(time_steps, n_results, 1, "Total boson number"), # here the number 1 refers to the number of trajectories.
    #     LineInfo(time_steps, n1_results, 1, "First site boson number") 
    # ]
    # save_data_path = joinpath(@__DIR__ , "plot_data")
    # makeplot(save_data_path, plot_lines...; xlabel = "t", ylabel = "I dont know", d, L, dt, t, k, w, U, J, state = "|01>")
end

function add_vline_gamma(p, dt)
    gamma = p / dt
    vline!([1 / gamma])
end

function check_rabi(U, J, dt)
    w = sqrt(U^2 + 2*(J)^2) / (2 * pi)
    p = w * 2 * dt
    println(p)
    add_vline_gamma(p, dt)
end

@time example()