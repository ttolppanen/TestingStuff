using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using Printf

function f()
    d = 2; L = 8; dt = 0.02; k = 6
    state = zeroone(d, L)
    H = bosehubbard(d, L)
    err = krylov_error_estimate(dt, k, state, H)
    @sprintf("err = %.3E", err)
end

f()