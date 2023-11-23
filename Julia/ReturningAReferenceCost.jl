using LinearAlgebra

function test_func(v)
    normalize!(v)
end

function f()
    v = rand(100)
    N = 1000
    @time for _ in 1:N
        test_func(v)
    end
end

f()