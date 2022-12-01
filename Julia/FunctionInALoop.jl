using LinearAlgebra

function calc(m)
    return m * m
end

function f()
    m = rand(3, 3)
    a = [1 2 3; 3 3 3; 4 4 4]
    test() = calc(a)
    m1 = test()
    @time for _ in 1:100
        m .= test()
    end
    @show norm(test() - m1)
    #@show m
end

f()