using LinearAlgebra
using SparseArrays

function matrix_vector(N)
    M = rand(N, N)
    v = rand(N)

    M[1:Int(end * 3 / 4), 1:Int(end * 3 / 4)] .= 0
    v[Int(end / 4):end] .= 0

    @time M * v
    return 
end

function sparse_matrix_vector(N)
    M = rand(N, N)
    v = rand(N)

    M[1:Int(end * 3 / 4), 1:Int(end * 3 / 4)] .= 0
    v[Int(end / 4):end] .= 0
    M = sparse(M)

    @time M * v
    return 
end

function matrix_sparse_vector(N)
    M = rand(N, N)
    v = rand(N)

    M[1:Int(end * 3 / 4), 1:Int(end * 3 / 4)] .= 0
    v[Int(end / 4):end] .= 0
    v = sparse(v)

    @time M * v
    return 
end

function sparse_matrix_sparse_vector(N)
    M = rand(N, N)
    v = rand(N)

    M[1:Int(end * 3 / 4), 1:Int(end * 3 / 4)] .= 0
    v[Int(end / 4):end] .= 0
    M = sparse(M)
    v = sparse(v)

    @time M * v
    return 
end

N = 10000
matrix_vector(N)
sparse_matrix_vector(N)
matrix_sparse_vector(N)
sparse_matrix_sparse_vector(N)

function matrix_structure_vector(N)
    M = rand(N, N)
    v = rand(N)

    v[1:Int(3 * N / 4)] .= 0
    M = sparse(M)
    v = sparse(v)

    @time for _ in 1:100
        M * v
    end
    return 
end

function structure_matrix_structure_vector(N)
    M = rand(N, N)
    v = rand(N)

    M[:, 1:Int(3 * N / 4)] .= 0
    v[1:Int(3 * N / 4)] .= 0
    M = sparse(M)
    v = sparse(v)
    
    @time for _ in 1:100
        M * v
    end
    return 
end

function cut_comparison(N)
    M = rand(Int(N / 4), Int(N / 4))
    v = rand(Int(N / 4))

    @time for _ in 1:100
        M * v
    end
    return 
end

function sparse_cut_comparison(N)
    M = rand(Int(N / 4), Int(N / 4))
    v = rand(Int(N / 4))

    M = sparse(M)
    v = sparse(v)

    @time for _ in 1:100
        M * v
    end
    return 
end

matrix_structure_vector(N)
structure_matrix_structure_vector(N)
cut_comparison(N)
sparse_cut_comparison(N)