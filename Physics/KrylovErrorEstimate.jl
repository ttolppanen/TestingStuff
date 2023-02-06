using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Printf

function f()
    d = 2; L = 8; dt = 0.32; k = 12
    state = zeroone(d, L)
    H = bosehubbard(d, L; w=0)
    err = krylov_error_estimate(dt, k, state, H)
    @sprintf("err = %.3E", err)
end

f()