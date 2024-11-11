println("")
include("EnvThreads.jl")
include("DisablePlots.jl")

println("")
println("Started!")
println("")

include("Physics/LeakageRemovalUnit/PlotsForPaper/second_figure/total_anh/gamma_rnd_msr/L2/code.jl")
include("Physics/LeakageRemovalUnit/PlotsForPaper/second_figure/total_anh/gamma_rnd_msr/normal_msr_comparison/code.jl")

println("")
println("Ended!")
