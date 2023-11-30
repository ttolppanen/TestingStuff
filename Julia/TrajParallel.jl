# function solvetrajectories(f::Function, traj::Integer; paral::Symbol = :none)
#     out = Array{Any}(undef, traj)
#     if paral == :threads
#         blas_threads = BLAS.get_num_threads()
#         BLAS.set_num_threads(1)
#         Threads.@threads for i in 1:traj
#             out[i] = f()
#         end
#         BLAS.set_num_threads(blas_threads)
#     elseif paral == :distributed
#         blas_threads = BLAS.get_num_threads()
#         BLAS.set_num_threads(1)
#         p_f(x) = f()
#         return pmap(p_f, 1:traj)
#         BLAS.set_num_threads(blas_threads)
#     elseif paral == :none
#         for i in 1:traj
#             out[i] = f()
#         end
#     else
#         throw(ArgumentError("possible keywords for paral {:none, :threads, :distributed}"))
#     end
#     return out
# end

function f()
    traj = 10000
    n_th = Threads.nthreads()
    trajs_to_calc = [traj ÷ n_th for _ in 1:n_th]
    trajs_to_calc[1] += traj % n_th
    m = rand(6, 6)
    
    out = [zeros(6,6) for _ in 1:traj]
    @sync for (i, traj_i) in enumerate(trajs_to_calc)
        Threads.@spawn begin
            trajs_calculated = 0
            if i != 1
                trajs_calculated = sum(trajs_to_calc[1:(i - 1)])
            end
            trajs_to_calculate = sum(trajs_to_calc[1:i])
            for t in (trajs_calculated + 1):trajs_to_calculate
                for _ in 1:100
                    exp(m)
                end
                out[t] .= exp(m)
            end
        end
    end
end

@time f()