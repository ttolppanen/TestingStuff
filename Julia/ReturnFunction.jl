function returnfunc(s)
    return asd(d) = d + s
end

function f()
    a = returnfunc(2)
    b = returnfunc(3)
    @show a(5)
end

f()