using SharedArrays

@everywhere function make_traj(mat)
    mat .= myid() * ones(2,2)
end

function f()
    N = 10
    a = SharedArray{Float64}(N)
    pre_v = []
    for i in 1:nprocs()
        push!(pre_v, @spawnat i zeros(N))
    end
    @sync @distributed for i = 1:N
        @spawnat myid() pre_v[myid()][i] = 1
    end

    display(fetch.(pre_v))
end

f()