include("./PAPER_ReadData.jl")
using MIPTDataProcessing
using QuantumStates
using QuantumOperators
using QuantumTimeEvolution

function f()
    # oldlayer = MakePlotFromOldData()
    filename = "mipt_same_results"
    layers = [getlayer("mipt_9", colorname = "w = 100"), getlayer("mipt_8", colorname = "w = 0")]
    prob = get_probabilities(filename, "mipt_9")
    title = get_attributes_string(filename, "mipt_9")
    p = plot(layers..., Guide.title(title))
    p |> SVG("test.svg", 6inch, 4inch)
end

function getlayer(groupname; colorname = "New Data")
    filename = "mipt_same_results"
    prob = get_probabilities(filename, groupname)
    dataend = get_data_end_mean(filename, groupname, prob, :entanglement_half)
    newLayer = layer(x=prob, y=dataend, Geom.line, color = [colorname])
end

function krylov(d, L, dt, t, k, traj, prob; color = :auto, kwargs...)
    function make_msrop(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op1 = n_bosons_projector(d, 1)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        return msrop
    end
    state0 = zeroone(d, L)
    H = bosehubbard(d, L; kwargs...)
    msrop = measurementoperators(nop(d), L)
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    calc_ent(r_traj) = trajmean(r_traj, s -> entanglement(d, L, s, Int(floor(L/2))))
    res = mipt(state0, H, k, meffect!, dt, t, prob, traj, calc_ent; paral = :distributed)
    return layer(x = prob, y = res, Geom.line, color = ["no data processing"])
end

f()