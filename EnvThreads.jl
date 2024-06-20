using Base.Threads
println("Threads used by julia = $(nthreads())")
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.resolve()