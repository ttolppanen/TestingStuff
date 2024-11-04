println("")
include("EnvThreads.jl")
include("DisablePlots.jl")

println("")
println("Started!")
println("")

include("Physics/LeakageRemovalUnit/optimal_J/equal_T1T2_LargeT2/T_times_L3/make_T_comparisons.jl")
include("Physics/LeakageRemovalUnit/optimal_J/equal_T1T2/T_times/make_T_comparisons.jl")
include("Physics/LeakageRemovalUnit/optimal_J/equal_T1T2/T_times_L3/make_T_comparisons.jl")

println("")
println("Ended!")
