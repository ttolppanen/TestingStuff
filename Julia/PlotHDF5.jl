using HDF5
using Plots

function same_param(g, params)
    if issetequal(keys(attributes(g)), String.(keys(params))) & all([read_attribute(g, key) == params[Symbol(key)] for key in keys(attributes(g))])
        return true
    end
    return false
end
function find_mipt_g(file; params...)
    for g_mipt in file
        if same_param(g_mipt, params)
            return g_mipt
        end
    end
    throw(ArgumentError("Couldn't find a group with those parameters"))
end
function f()
    x = 0:0.01:0.4
    y =[]
    h5open("mipt_test_data.h5", "r") do file
        g_mipt = find_mipt_g(file; d=2, L=6, dt=0.02, t=30, k=6, initial_state = "zeroone(d, L)", w=1, U=1, J=1, msr_type = "std_n")
        for p in x
            push!(y, read_attribute(g_mipt["p = $p/entanglement_h"], "end_mean"))
        end
    end
    plot(x, y) 
end

@time f()