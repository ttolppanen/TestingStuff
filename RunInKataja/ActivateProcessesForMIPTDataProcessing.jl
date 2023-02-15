using Distributed
@everywhere begin
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
end
@everywhere using MIPTDataProcessing