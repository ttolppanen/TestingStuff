function f()
    dt = 0.02
    d = 4
    prob = 0.15:0.015:0.9
    traj = 50
    observables = 
    [
        :entanglement_half,
        :half_n_fluctuations
    ]
    for maxdim in [20, 40]
        for L in 8:2:12
            @time generate_data("mipt_same_results"; alg= :mps, d, L, dt, t=30, state= :zeroone, trotter_order=2, 
            probabilities=prob, measurement= :std_n, traj, w=0, U=5, save_before_effect = true, 
            J=1, observables, cutoff = 1E-8, maxdim)
        end
    end
end

f()