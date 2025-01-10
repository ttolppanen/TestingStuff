using PlotAndSave

function f()
    path1 = joinpath(@__DIR__, "gamma")
    path2 = joinpath(@__DIR__, "gamma_tail")
    new_path = joinpath(@__DIR__, "gamma_with_tail")
    merge_lines(new_path, path1, path2)
end

f()