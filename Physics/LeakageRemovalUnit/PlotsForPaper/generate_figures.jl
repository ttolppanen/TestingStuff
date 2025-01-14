# Figure2
include(joinpath(@__DIR__, "2_figure/first_panel/time/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "2_figure/second_panel/time/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "2_figure/big_panel/combined/paper_plots/paper_plot.jl"))

# Figure3
include(joinpath(@__DIR__, "3_figure/disorder/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "3_figure/disorder_realistic/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "3_figure/J/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "3_figure/J_realistic/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "3_figure/size/paper_plots/paper_plot.jl"))
include(joinpath(@__DIR__, "3_figure/size_realistic/paper_plots/paper_plot.jl"))