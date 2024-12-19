using QuantumOperators
using QuantumStates
using PlotAndSave

function f()
    d = 3
    L = 6 - 1
    param = Dict()
    param[:J] = ValueInfo(5, "5MHz")
    param[:w] = ValueInfo(7500 / param[:J].val, "7.5GHz")
    param[:dw] = ValueInfo(100 / param[:J].val, "100MHz")
    param[:U] = ValueInfo(250 / param[:J].val, "250MHz")
    param[:T] = ValueInfo(0.1, "0.1K")

    t_T = param[:T].val
    t_w = param[:w].val * param[:J].val * 2 * pi * 10^6
    t_U = param[:U].val * param[:J].val * 2 * pi * 10^6
    
    n_op = nall(d, L)
    N = 1000
    out = []
    for _ in 1:N
        ts_sample = sample_transmon_thermal_state(d, L, t_T, t_w, t_U)
        push!(out, expval(ts_sample, n_op))
    end
    sum(out) / N
end
f()