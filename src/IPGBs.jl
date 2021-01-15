#TODO re-export the relevant functions: the 4ti2 interface and my GB
#implementations
module IPGBs

include("./FastBitSets.jl")
include("./GBElements.jl")
include("./SupportTrees.jl")
include("./Binomials.jl")
include("./GradedBinomials.jl")
include("./GBTools.jl")

include("./SignaturePolynomials.jl")
include("./CriticalPairs.jl")
include("./BinomialSets.jl")
include("./GBAlgorithms.jl")

include("./Buchberger.jl")
include("./SignatureAlgorithms.jl")
include("./FourTi2.jl")

using GBAlgorithms
using Buchberger
using SignatureAlgorithms

#TODO Finish implementing these. They should be practical ways to call my
#algorithms
function buchberger_algorithm(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    algorithm = BuchbergerAlgorithm()
    run(algorithm)
    return current_basis(algorithm)
end

#TODO this looks almost like the one above but duplicated. Maybe there's a better
#way to do this
function signature_algorithm(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    algorithm = SignatureAlgorithm()
    run(algorithm)
    return current_basis(algorithm)
end

end
