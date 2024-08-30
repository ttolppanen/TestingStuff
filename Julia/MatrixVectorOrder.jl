function f()
    n = 100
    m = rand(n, n)
    v = rand(n)

    @time (m * m) * v
    @time m * m * v
    @time m * (m * v)
    @time begin
        s = m * v
        m * s
    end
end

f();