using QuantumStates
using ITensors

function f()
    @show Threads.nthreads()
    out = [[] for _ in 1:Threads.nthreads()]
    m = rand(4,4)
    Threads.@threads for i in 1:30
        # @show Threads.threadid()
        s = Threads.threadid() + 1
        @show i
        push!(out[Threads.threadid()], deepcopy(m))
    end
    @show sum(reduce(vcat, out))
end

function mps_thread()
    d = 4; L = 5
    state = zeroonemps(d, L)
    nop = op("N", siteinds(state), 1)
    Threads.@threads for _ in 1:30
        s = deepcopy(state)
        s = apply(nop, state)
        normalize!(s)
    end
    state
end

mps_thread()