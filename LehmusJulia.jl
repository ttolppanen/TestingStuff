println("")
include("EnvThreads.jl")
include("DisablePlots.jl")

println("")
println("Started!")
println("")

include("Physics/LeakageRemovalUnit/optimal_J/T_times/make_T_comparisons.jl")

println("")
println("Ended!")
