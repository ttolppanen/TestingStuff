println("\n Updating environment...\n")
using Pkg
Pkg.Registry.update()
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
Pkg.precompile()
println("\n Done! \n")