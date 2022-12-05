using LinearAlgebra

function f()
    d = 2^12
    v = rand(d)
    m = rand(d, d) * 0.0001
    out = [deepcopy(v)]
    @time for _ in 1:10
        push!(out, deepcopy(out[end]))
        mul!(out[end], m, out[end])
    end
    k = 100
    @time out = [deepcopy(v) for _ in 1:k]
    # @time out = [zeros(size(v)[1]) for _ in 1:k]
    out[1][1] = 69
    @show out[2][1]
end

f();