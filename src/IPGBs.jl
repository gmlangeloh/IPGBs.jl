#TODO re-export the relevant functions: the 4ti2 interface and my GB
#implementations
module IPGBs
export groebner_basis

include("./SolverTools.jl")
include("./FastBitSets.jl")
include("./IPInstances.jl")
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
include("./Markov.jl")

using .GBAlgorithms
using .IPInstances
import .Buchberger: BuchbergerAlgorithm
import .SignatureAlgorithms: SignatureAlgorithm
import .Binomials: Binomial
import .GradedBinomials: GradedBinomial

"""
Simple truncation works when the data in A, b are all non-negative.
This is more efficient than IP truncation, and is exact, contrary to LP
truncation, so it should be used whenever possible.

This function returns whether simple truncation can be used.
"""
function use_simple_truncation(
    A :: Array{Int, 2},
    b :: Vector{Int}
) :: Bool
    m, n = size(A)
    return all(A[i, j] >= 0 for j in 1:n for i in 1:m) &&
        all(b[i] >= 0 for i in 1:m)
end

function parse_truncation(
    truncation_type :: Symbol,
    A :: Array{Int, 2},
    b :: Vector{Int}
) :: Tuple{Symbol, DataType}
    var_type = Any #Only matters when some LP/IP model will be built
    if truncation_type == :Heuristic
        #Choose either LP or simple, depending on the form of A and b
        if use_simple_truncation(A, b)
            truncation_type = :Simple
        else
            truncation_type = :LP
        end
    end
    if truncation_type == :LP
        type = :Model
        var_type = Real
    elseif truncation_type == :IP
        type = :Model
        var_type = Int
    elseif truncation_type == :Simple
        type = :Simple
    elseif truncation_type == :None
        type = :None
        var_type = Any
    else
        error("Unknown truncation type.")
    end
    return type, var_type
end

"""
TODO write full documentation for this top-level function...

TODO as the number of arguments and complexity of preprocessing goes up, that
should probably be done in another module

- module_order is one of :ltpot, :pot and :top. Only relevant for signature
algorithms
- truncation type is one of :None, :Heuristic, :LP, :IP, :Simple
"""
function groebner_basis(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{T, 2},
    u :: Vector{Int};
    use_signatures :: Bool = false,
    implicit_representation :: Bool = false,
    module_order :: Symbol = :ltpot,
    truncation_type :: Symbol = :Heuristic,
    quiet :: Bool = false,
    minimization :: Bool = true
) :: Vector{Vector{Int}} where {T <: Real}
    #Setting parameters
    algorithm_type = use_signatures ? SignatureAlgorithm : BuchbergerAlgorithm
    representation = implicit_representation ? GradedBinomial : Binomial
    use_minimization = implicit_representation ? false : minimization
    trunc_type, trunc_var = parse_truncation(truncation_type, A, b)
    instance = IPInstance(A, b, C, u)

    #Signatures aren't currently supported with implicit representation
    @assert !(use_signatures && implicit_representation)

    #Run GB algorithm over the given instance
    if use_signatures
        algorithm = algorithm_type(
            representation, obj, A, b, module_order, truncate, use_minimization
        )
    else
        algorithm = algorithm_type(
            representation, instance, trunc_type, trunc_var, use_minimization
        )
    end
    results = GBAlgorithms.run(algorithm, instance, quiet=quiet)
    return results
end

end
