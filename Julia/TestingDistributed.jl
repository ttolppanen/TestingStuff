using Distributed
if length(procs()) == 1
    addprocs()
end
@show length(procs())

@everywhere function run_f(x)
    rand(10000)
end 

function f()
    traj = 10000
    @time res = pmap(run_f, 1:traj)
    #@show res
end

f();