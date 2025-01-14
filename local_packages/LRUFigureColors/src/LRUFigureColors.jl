module LRUFigureColors

using Plots

export color_small_msr
export color_high_msr
export color_periodic
export color_dissipation_removal
export color_t2
export color_t1
export color_ts

color_small_msr = RGB(245 /255, 76 /255, 181 /255)
color_high_msr = RGB(91 /255, 245/255, 76/255)
color_periodic = RGB(76 /255, 146 /255, 245 /255)
color_dissipation_removal = RGB(245 /255, 187 /255, 76 /255)
color_t2 = RGB(224 /255, 154 /255, 22 /255)
color_t1 = RGB(224 /255, 203 /255, 22 /255)
color_ts = color_periodic

# color_small_msr = :green
# color_high_msr = :red
# color_periodic = :blue
# color_dissipation_removal = :yellow
# color_t2 = :orange
# color_t1 = :yellow
# color_ts = color_periodic

end # module FigureColors
