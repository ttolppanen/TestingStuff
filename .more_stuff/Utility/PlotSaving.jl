using Plots

function filename()
    name = split(Base.source_path(), "\\")[end]
    name = chop(name, head=0, tail=3) #remove .jl
    return name
end

function saveplot(pl, name)
    path = joinpath(pwd(), "Physics", "plots", filename())
    mkpath(path)
    savefig(pl, joinpath(path, name * ".png"))
end
@show joinpath(dirname(@__FILE__), "pltohere.png")