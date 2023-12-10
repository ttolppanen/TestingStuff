using ExponentialUtilities
using QuantumStates
using QuantumOperators
using QuantumTimeEvolution
using SparseArrays
using StaticArrays
using LinearAlgebra

function f()
    d = 2; L = 10;
    dt = 0.02
    k = 6
    state = Vector(zeroone(d, L))
    H = bosehubbard(d, L)
    cache = ExpvCache{Float64}(k)
    Ks = KrylovSubspace{ComplexF64, Float64}(d^L, k)
    pa_k = PA_krylov(d^L, k)
    take_time_step! = take_krylov_time_step_function(H, k, dt, pa_k)
    take_time_step_new! = (state) -> begin
        lanczos!(Ks, H, state)
        expv!(state, dt, Ks; cache)
        normalize!(state)
    end

    test_mine(take_time_step!, state)
    @time test_mine(take_time_step!, state)

    test_expkit(Ks, H, state, dt, cache)
    @time test_expkit(Ks, H, state, dt, cache)

    test_mine(take_time_step_new!, state)
    @time test_mine(take_time_step_new!, state)
end

function test_expkit(Ks, H, state, dt, cache)
    for _ in 1:1500
        lanczos!(Ks, H, state)
        expv!(state, dt, Ks; cache)
        normalize!(state)
    end
end

function test_mine(f, state)
    for _ in 1:1500
        state .= f(state)
    end
end

function take_krylov_time_step_function(H::AbstractMatrix{<:Number}, k::Integer, dt::Real, pa_k::PA_krylov)
    let
        H = H
        k = k
        dt = dt
        pa_k = pa_k
        function take_time_step!(state)
            QuantumTimeEvolution.krylovsubspace!(state, H, k, pa_k) # makes changes into pa_k
            if !all(isfinite, pa_k.H_k) 
                throw(ArgumentError("Hₖ contains Infs or NaNs. This is is usually because k is too small, too large or there is no time evolution H * state0 = 0.")) 
            end
            mul!(state, pa_k.U, @view(exp(-1im * dt * pa_k.H_k)[:, 1]))
            normalize!(state)
            return state
        end
        return take_time_step!
    end
end

f()