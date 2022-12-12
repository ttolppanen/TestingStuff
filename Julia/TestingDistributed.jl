using Distributed
addprocs(8)

@everywhere function run_f(x)
    rand(10000)
end 

function f()
    traj = 10
    @time res = pmap(run_f, 1:traj)
    #@show res
end

f();