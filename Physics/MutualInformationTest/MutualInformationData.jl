using MIPTDataProcessing

function make_data()
    # dts = [0.02, 0.04, 0.08]
    # probs = [0.01:0.005:0.06, 0.02:0.01:0.12, 0.04:0.02:0.24]
    # traj = 3000
    # @time generate_data("mipt_same_results"; alg= :krylov, d=2, L=8, dt=dts[1], t=30, state= :zeroone, k=6, 
    # probabilities=probs[1], measurement= :std_n, traj, w=0, U=10, 
    # J=1, observables=[:entanglement_half])
    # @time generate_data("mipt_same_results"; alg= :krylov, d=2, L=8, dt=dts[2], t=30, state= :onezero, k=6, 
    # probabilities=probs[2], measurement= :std_n, traj, w=0, U=10, 
    # J=1, observables=[:entanglement_half])
    # @time generate_data("mipt_same_results"; alg= :krylov, d=2, L=8, dt=dts[3], t=30, state= :onezero, k=6, 
    # probabilities=probs[3], measurement= :std_n, traj, w=0, U=10, 
    # J=1, observables=[:entanglement_half])
end

function make_collapse_data()
    dt = 0.02
    prob = 0.001:0.001:0.02
    traj = 5000
    for L in 6:2:10
        println("first")
        @time generate_data("mipt_same_results"; alg= :krylov, d=2, L, dt, t=30, state= :zeroone, k=6, 
        probabilities=prob, measurement= :std_n, traj, w=0, U=10, save_before_effect = false, 
        J=1, observables=[:entanglement_half, :mutual_information_first_to_half])
        println("second")
        @time generate_data("mipt_same_results"; alg= :krylov, d=2, L, dt, t=30, state= :zeroone, k=6, 
        probabilities=prob, measurement= :std_n, traj, w=0, U=10, save_before_effect = true, 
        J=1, observables=[:entanglement_half, :mutual_information_first_to_half])
    end
end

# make_data()
# make_data_higher_k()
make_collapse_data()