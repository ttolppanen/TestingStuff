using Distributed
addprocs(16)

@everywhere function run_f(x)
    r_f() = rand(10000)
    r_f()
end 

function f()
    traj = 1000
    a = Array{Any}(undef, traj)
    @time for i in 1:traj
        a[i] = r_f()
    end
    @time res = pmap(run_f, 1:traj)
    #@show res
end

f();
f();