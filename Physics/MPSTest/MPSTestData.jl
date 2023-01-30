using MIPTDataProcessing

function make_data()
    dt = 0.3
    prob = 0.15:0.075:0.9
    traj = 500
    for L in 8:2:14
        @time generate_data("mipt_same_results"; alg= :mps, d=2, L, dt=0.3, t=20, state= :zeroone, trotter_order=2, 
        probabilities=prob, measurement= :std_n, traj, w=0, U=10, maxdim = 10, cutoff=1E-8,
        J=1, observables=[:entanglement_half])
    end
end

make_data()