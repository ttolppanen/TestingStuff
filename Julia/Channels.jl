using LinearAlgebra

function f()
    traj = 100
    size = 1000
    iterations = 30
    m1 = zeros(size, size)
    m2 = zeros(size, size)
    out = Channel{Matrix{Float64}}(traj)
    out2 = Channel{Matrix{Float64}}(traj)
    for i in 1:traj
        put!(out, i * ones(size, size))
    end
    for _ in 1:traj
        m = take!(out)
        # display(m);
        put!(out, m)
    end
    @time @sync for i in 1:traj
        Threads.@spawn begin
            m = take!(out)
            for _ in 1:iterations
                mul!(m, m1, m2)
            end
            put!(out2, m)
        end
    end
    # for _ in 1:traj
    #     m = take!(out)
    #     # display(m);
    #     put!(out, m)
    # end
end


f()