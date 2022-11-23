using QuantumStates
using QuantumOperators
using LinearAlgebra

#entanglement with schmidt is significantly faster

function f()
    d = 4; L = 8
    state = normalize(zeroone(d,L) + onezero(d, L))
    res = @time entanglement(d, L, state, Int(L/2))
    res2 = @time entanglement_sch(d, L, state, Int(L/2))
    @show res - res2
end

@profview f()