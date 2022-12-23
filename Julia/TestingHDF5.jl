using HDF5
using Dates
using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Statistics

function same_param(g, params)
    if issetequal(keys(attributes(g)), String.(keys(params))) & all([read_attribute(g, key) == params[Symbol(key)] for key in keys(attributes(g))])
        return true
    end
    return false
end
function saveh5(ent, msr_prob; simulation_param...)
    h5open("mipt_test_data.h5", "cw") do file
        for g_mipt in file
            if same_param(g_mipt, simulation_param)
                if haskey(g_mipt, "p = $msr_prob")
                    println("Data already exists for p = $msr_prob")
                    return
                else
                    g_p = create_group(g_mipt, "p = $(msr_prob)")
                    attributes(g_p)["p"] = msr_prob   
                    data_size = size(ent) # 1 timesteps, 2 trajectories
                    data = create_dataset(g_p, "entanglement_h", Float64, ((data_size[1], data_size[2]), (data_size[1], -1)), chunk = (data_size[1], 1))
                    data[:, :] = ent
                    attributes(data)["end_mean"] = mean(ent[end, :])
                end
                return
            end
        end
        # make new

        group_name = "mipt_" * string(length(file) + 1)
        g_mipt = create_group(file, group_name)
        for (key, value) in simulation_param
            attributes(g_mipt)["$key"] = value
        end
        g_p = create_group(g_mipt, "p = $(msr_prob)")
        attributes(g_p)["p"] = msr_prob   
        data_size = size(ent) # 1 timesteps, 2 trajectories
        data = create_dataset(g_p, "entanglement_h", Float64, ((data_size[1], data_size[2]), (data_size[1], -1)), chunk = (data_size[1], 1))
        data[:, :] = ent
        attributes(data)["end_mean"] = mean(ent[end, :])
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
    saveh5(ent_data, msr_prob; d, L, dt, t, k, initial_state = "zeroone(d, L)", w=1, U=1, J=1, msr_type = "std_n")
end

function read_data()
    h5open("mipt_test_data.h5", "r") do file
        ent = file["mipt_22_12_2022/entanglement_h"][end, :]
        @show sum(ent) / 100
    end
end

begin
    for L in [4, 6, 8]
        for p in 0:0.001:0.1
            @time calc_ent(2, L, 0.02, 30, 6, 100, p);
        end
    end
end