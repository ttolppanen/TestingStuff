using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere Pkg.instantiate()