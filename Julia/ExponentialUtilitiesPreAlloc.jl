using ExponentialUtilities
using LinearAlgebra

function f()
    d = 1000; t = -1im * 0.1; k = 6; # krylov dimension
    v = complex(normalize(rand(d)))
    m = rand(d, d) + 1.0im * rand(d, d)
    m = m + m'
    println(ishermitian(m))

    cache = ExpvCache{Float64}(k)
    ks = KrylovSubspace{ComplexF64, Float64}(d, k)

    test_exp_utl(ks, m , v, t; cache)
    v = complex(normalize(rand(d)))

    @time test_exp_utl(ks, m, v, t; cache)
    v = complex(normalize(rand(d)))
    @time test_exp_utl(ks, m, v, t; cache)
    v = complex(normalize(rand(d)))
    @time test_exp_utl(ks, m, v, t; cache)
end

function test_exp_utl(ks, m, vec, t; kwargs...) # kwargs for cache
    lanczos!(ks, m, vec)
    expv!(vec, t, ks; kwargs...)
    normalize!(vec)
end

f();