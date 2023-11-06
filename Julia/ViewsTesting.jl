using LinearAlgebra

function change_matrix(m)
    s = size(m)
    for i in 1:s[1]
        for j in 1:s[2]
            m[i, j] = 1
        end
    end
end

function f()
    m = rand(5, 5)
    @time change_matrix(@view m[2:4, 2:4])
    m
end

f()