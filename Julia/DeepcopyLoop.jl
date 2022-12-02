using LinearAlgebra

function f()
    d = 100
    v = rand(d)
    m = rand(d, d) * 0.0001
    out = [deepcopy(v)]
    @time for _ in 1:100
        push!(out, deepcopy(out[end]))
        mul!(out[end], m, out[end])
    end
end

f()