function f()
    @show Threads.nthreads()
    out = [[] for _ in 1:Threads.nthreads()]
    m = rand(4,4)
    Threads.@threads for _ in 1:30
        @show Threads.threadid()
        push!(out[Threads.threadid()], deepcopy(m))
    end
    @show sum(reduce(vcat, out))
end

f()