#TODO re-export the relevant functions: the 4ti2 interface and my GB
#implementations
module IPGBs
export groebner_basis

include("./FastBitSets.jl")
include("./Statistics.jl")
include("./GBTools.jl")
include("./Orders.jl")
include("./FastComparator.jl")
include("./GBElements.jl")
include("./SupportTrees.jl")
include("./Binomials.jl")
include("./GradedBinomials.jl")

include("./BinomialSets.jl")
include("./SignaturePolynomials.jl")
include("./TriangleHeaps.jl")
include("./GBAlgorithms.jl")

include("./Buchberger.jl")
include("./SignatureAlgorithms.jl")
include("./FourTi2.jl")

using .GBAlgorithms
import .Buchberger: BuchbergerAlgorithm
import .SignatureAlgorithms: SignatureAlgorithm
import .Binomials: Binomial
import .GradedBinomials: GradedBinomial

function groebner_basis(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{T, 2},
    u :: Vector{Int};
    use_signatures :: Bool = false,
    implicit_representation :: Bool = false,
    module_order :: Symbol = :ltpot,
    truncate :: Bool = true,
    quiet :: Bool = false,
    minimization :: Bool = true
) :: Vector{Vector{Int}} where {T <: Real}
    #Setting parameters
    algorithm_type = use_signatures ? SignatureAlgorithm : BuchbergerAlgorithm
    representation = implicit_representation ? GradedBinomial : Binomial
    use_minimization = implicit_representation ? false : minimization
    obj = Float64.(C) #Consider the objective function as floating point regardless

    #Signatures aren't currently supported with implicit representation
    @assert !(use_signatures && implicit_representation)

    #Run GB algorithm over the given instance
    if use_signatures
        algorithm = algorithm_type(representation, obj, A, b, module_order, truncate, use_minimization)
    else
        algorithm = algorithm_type(representation, obj, A, b, truncate, use_minimization)
    end
    results = GBAlgorithms.run(algorithm, A, b, obj, u, quiet=quiet)
    return results
end

end
