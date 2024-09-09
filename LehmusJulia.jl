println("")
include("EnvThreads.jl")
include("DisablePlots.jl")

println("")
println("Started!")
println("")

include("Physics/LeakageRemovalUnit/PlotsForPaper/fourth_figure/disorder/make_T_comparisons.jl")

println("")
println("Ended!")
