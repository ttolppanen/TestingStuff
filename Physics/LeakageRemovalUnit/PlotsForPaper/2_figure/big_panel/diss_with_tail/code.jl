using PlotAndSave

function f()
    path1 = joinpath(@__DIR__, "diss")
    path2 = joinpath(@__DIR__, "diss_tail")
    new_path = joinpath(@__DIR__, "diss_with_tail")
    merge_lines(new_path, path1, path2)
end

f()