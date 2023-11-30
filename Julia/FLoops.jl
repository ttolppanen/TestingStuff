using FLoops

function f()
    N = 10
    out = [zeros(2,2) for _ in 1:N]
    @floop for i in 1:N
        out[i] .= i * v
    end
    display(out)
end

f()