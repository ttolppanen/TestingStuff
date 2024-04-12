using Distributed
using JLD2
using Plots
using LinearAlgebra
using SparseArrays
@everywhere using QuantumStates
@everywhere using QuantumOperators
using PlotAndSave
@everywhere using QuantumTimeEvolution
using OrderedCollections
using DSP

function traj_mean(result)
    out = zeros(size(result[1]))
    for traj in result
        out .+= traj
    end
    return out ./ length(result)
end

function qubit_error_fix_U(d, L, param, dw, dU)
    H = sparse(complex(zeros(9, 9)))
    H[d + 1, 2] = param[:J]
    H[2, d + 1] = param[:J]

    U = sparse(exp(1im * Matrix(H) * param[:dt]))
    sub_d = d^(param[:L]- L)
    return kron(U, sparse(I, sub_d, sub_d))
end

function measure_ends_effect!(msr_op::Vector, 
    dt, msr_rate, feedback::Vector)
    time = [rand() / msr_rate]
    out(state) = begin
        time[1] += dt

        if time[1] < 1 / msr_rate
            return state
        end
        time[1] -= 1 / msr_rate

        for site in [length(msr_op)]
            measuresite!(state, msr_op, site)
            state .= feedback[site] * state
            normalize!(state)
        end
        return state
    end
    return out
end

function initial_parameters(L)
    param = OrderedDict()
    param[:gamma] = 27
    param[:dt] = 1
    param[:t] = 400
    param[:d] = 3
    param[:L] = L
    state = normalize(bosonstack(1, L, 1; d = param[:d]) + allzero(param[:d], L))
    param[:state] = ValueInfo(state, "|10..> + |00..>")
    param[:alg] = "exact"
    param[:H] = "BH, d_w, fixed |2> energy"
    param[:extra_U] = "two qubit sigma_x"
    param[:w] = 1000
    param[:dw] = 10
    param[:U] = 50
    param[:J] = 1

    d = param[:d]
    msrop = measurementoperators(nop(d), L)
    feedback = [1 1 1; 0 0 0; 0 0 0]
    fb_op = [singlesite(feedback, L, i) for i in 1:L]
    
    param[:effect] = ValueInfo((msr = msrop, fb = fb_op), "rand start, msr_n last site -> fb -> 0")
    return param
end

function f()
    param = initial_parameters(4)
    # -----------
    traj = 200
    Ls = [3,4,5]
    # -----------
    lines = []
    for L in Ls
        param = initial_parameters(L)
        d = param[:d]
        state = param[:state].val

        msrop = param[:effect].val.msr
        feedback = param[:effect].val.fb

        sigma_x = sparse(complex([0 1.0 0; 1 0 0; 0 0 0]))
        ops = singlesite(sigma_x, L, 1)
        obs = (state) -> expval(state, ops)

        dt = param[:dt]; t = param[:t]
        pa_args, pa_out = pre_alloc_threads(zeros(1, length(0:dt:t)), traj, deepcopy(state))

        out = []

        r_f(out, work_vector) = begin
            dw = param[:w] .+ (rand(L) .- 0.5) * param[:dw]
            dU = 2 * dw .- (2 * param[:w] - param[:U])
            H = bosehubbard(d, L; w = dw, U = dU, J = param[:J])
            Uop = exp(-1im * dt * Matrix(H))

            error_U = qubit_error_fix_U(d, 2, param, dw, dU)
            fix_error(state) = begin
                state .= error_U * state
                return state
            end

            msr_effect! = measure_ends_effect!(msrop, param[:dt], param[:gamma], feedback)
            effect!(state) = begin
                state = fix_error(state)
                return msr_effect!(state)
            end
            exactevolve(state, work_vector, Uop, dt, t, obs; effect!, out)
        end
        solvetrajectories(r_f, traj, pa_args, pa_out; paral = :threads)
        r = map(x -> abs.(hilbert(x[1,:])), pa_out)
        r = traj_mean(r)

        push!(lines, LineInfo(0:dt:t, r, traj, "L = $L"))
    end
    path = joinpath(@__DIR__)
    makeplot(path, lines...; xlabel = "t", ylabel = "<first site sigma_x>", copy_code = false, param...)
    # addtrajectories(path, line)
end

@time f()