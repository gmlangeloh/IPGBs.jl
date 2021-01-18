#TODO re-export the relevant functions: the 4ti2 interface and my GB
#implementations
module IPGBs
export groebner_basis

include("./FastBitSets.jl")
include("./GBElements.jl")
include("./SupportTrees.jl")
include("./Binomials.jl")
include("./GradedBinomials.jl")
include("./GBTools.jl")

include("./BinomialSets.jl")
include("./SignaturePolynomials.jl")
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
    C :: Array{Int, 2},
    u :: Vector{Int};
    use_signatures :: Bool = false,
    implicit_representation :: Bool = false
) :: Vector{Vector{Int}}
    #Setting parameters
    algorithm_type = use_signatures ? SignatureAlgorithm : BuchbergerAlgorithm
    representation = implicit_representation ? GradedBinomial : Binomial

    #Signatures aren't currently supported with implicit representation
    @assert !(use_signatures && implicit_representation)

    #Run GB algorithm over the given instance
    algorithm = algorithm_type(representation, C)
    return GBAlgorithms.run(algorithm, A, b, C, u)
end

end
