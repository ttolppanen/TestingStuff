using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere Pkg.Registry.update()
@everywhere Pkg.resolve()
@everywhere Pkg.instantiate()