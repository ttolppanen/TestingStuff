using HDF5
using Dates
using QuantumStates
using QuantumOperators
using QuantumTimeEvolution

function saveh5(ent)
    h5open("mipt_test_data.h5", "w") do file
        group_name = "mipt_" * Dates.format(Dates.now(), "d_m_Y")
        g = create_group(file, group_name)
        data_size = size(ent) # 1 timesteps, 2 trajectories
        data = create_dataset(g, "entanglement_h", Float64, ((data_size[1], data_size[2]), (data_size[1], -1)), chunk = (data_size[1], 1))
        data[:, :] = ent
    end
end

function calc_ent(d, L, dt, t, k, traj, msr_prob)
    state = zeroone(d, L)
    H = bosehubbard(d, L)
    effect! = measuresitesrandomly!(L, nop(d), msr_prob)
    r_f() = krylovevolve(state, H, dt, t, k; effect!)
    r = solvetrajectories(r_f, traj)
    ent_data = zeros(length(r[1]), length(r))
    for traj_i in eachindex(r) 
        for timestep_i in eachindex(r[traj_i])
            s = r[traj_i][timestep_i]
            ent_data[timestep_i, traj_i] = entanglement(d, L, s, Int(floor(L/2)))
        end
    end
    @time saveh5(ent_data)
end

function read_data()
    h5open("mipt_test_data.h5", "r") do file
        ent = file["mipt_22_12_2022/entanglement_h"][end, :]
        @show sum(ent) / 100
    end
end

# calc_ent(2, 6, 0.02, 30, 6, 100, 0.1);
@time read_data()