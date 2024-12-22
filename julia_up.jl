println("\n Updating environment...\n")
using Pkg
Pkg.Registry.update()
Pkg.activate(".")
Pkg.instantiate()
Pkg.resolve()
Pkg.precompile()
println("\n Done! \n")