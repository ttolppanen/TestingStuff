function confirm_registry_update()
    println("Have you ran Pkg.Registry.update()? y/n")
    out = readline()
    if out != "y"
        println("Please restart julia just with a single process and run \n Using Pkg \n Pkg.Registry.update() \n and then restart julia with processes")
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