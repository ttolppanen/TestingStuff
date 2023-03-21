function confirm_registry_update()
    println("Have you ran Pkg.Registry.update()? y/n")
    out = readline()
    if out != "y"
        exit()
    end
end

confirm_registry_update()

using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere Pkg.resolve()
@everywhere Pkg.instantiate()