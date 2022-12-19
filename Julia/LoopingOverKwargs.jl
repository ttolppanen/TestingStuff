function tf(; kwargs...)
    out = ""
    for (key, value) in kwargs
        out *= "$key = $value, "
    end
    out[1:end-2]
end

begin
    dt = 0.2
    x = 55
    L = 7
    prob = 0:01:10
    tf(;dt, x, L, prob) 
end