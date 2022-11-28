struct test
    a
    b
end

function f()
    m1 = [1 2; 3 4]
    m2 = [2 2; 3 3]
    t = test(m1, m2)
    t.a[1,1] = 7
    @show t.a
    @show m1
end

f()