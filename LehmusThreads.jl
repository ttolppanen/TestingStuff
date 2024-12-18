using Base.Threads
println("Threads used by julia = $(nthreads())")

using Pkg # Activates environment
Pkg.activate(".")

using Plots # Disables plots
ENV["GKSwstype"] = "nul"

println("")
println("Started!")
println("")