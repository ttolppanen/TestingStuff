using LinearAlgebra

function test_func(a, b)
     return a, b
end

function test_this(a, b, N)
    t = (a, b)
    @code_warntype for _ in 1:N
        t = test_func(t...)
    end
end

function f()
    N = 10000
    a = rand(N)
    b = rand()
    @code_warntype test_this(a, b, N)
end

f()