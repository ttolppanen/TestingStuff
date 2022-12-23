using HDF5
using Statistics

function make_file(n_res, n_traj, prob)
    h5open("prob_test.h5", "w") do file
        g = create_group(file, "mipt")
        for p in prob
            g_p = create_group(g, "p = $p")
            attributes(g_p)["p"] = p
            data = create_dataset(g_p, "ent", Float64, ((n_res, n_traj), (n_res, -1)), chunk = (n_res, 1))
            data[:, :] = rand(n_res, n_traj)
        end
    end
end

function loop_prob()
    h5open("prob_test.h5", "r") do file
        for obj in file["mipt"]
            @show read_attribute(obj["ent"], "steady_state_mean")
        end
    end
end

# @time make_file(2500, 10000, 1:100)
@time loop_prob()