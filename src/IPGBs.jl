module IPGBs
export groebner_basis, markov_basis, optimize, gb_heuristic, gb_heuristic!, moip_gb_solve

include("./Globals.jl")
include("./GBTools.jl")
include("./MatrixTools.jl")
include("./SolverTools.jl")
include("./IPInstances.jl")
include("./LatticePrograms.jl")
include("./CombinatorialOptimizationInstances.jl")

include("./FastBitSets.jl")
include("./Statistics.jl")
include("./Orders.jl")
include("./MonomialHeaps.jl")
include("./FastComparator.jl")
include("./Pairs.jl")
include("./GBElements.jl")
include("./SupportTrees.jl")
include("./SupportMatrixTrees.jl")
include("./Binomials.jl")
include("./GradedBinomials.jl")

include("./BinomialSets.jl")
include("./SignaturePolynomials.jl")
include("./TriangleHeaps.jl")
include("./FourTi2.jl")
include("./GBAlgorithms.jl")

include("./ZDDs.jl")
include("./FGLM.jl")
include("./Buchberger.jl")
include("./SignatureAlgorithms.jl")
include("./Markov.jl")
include("./StandardDecomposition.jl")
include("./Optimize.jl")
include("./Walkback.jl")
include("./FeasibleGraphs.jl")

include("./SingleObjective.jl")
include("./MultiObjectiveTools.jl")
include("./MultiObjectiveStats.jl")
include("./MultiObjectiveAlgorithms.jl")

using .GBAlgorithms
using .GBTools
using .IPInstances
import .Markov: markov_basis
import .Optimize: optimize, optimize_with!, gb_heuristic, gb_heuristic!
import .Buchberger: BuchbergerAlgorithm
import .SignatureAlgorithms: SignatureAlgorithm
import .Binomials: Binomial
import .GradedBinomials: GradedBinomial
import .MultiObjectiveAlgorithms: moip_gb_solve

using GLPK
using JuMP
using Random

"""
Simple truncation works when the data in A, b are all non-negative.
This is more efficient than IP truncation, and is exact, contrary to LP
truncation, so it should be used whenever possible.

This function returns whether simple truncation can be used.
"""
function use_simple_truncation(
    A::Array{Int,2},
    b::Vector{Int}
)::Bool
    m, n = size(A)
    return all(A[i, j] >= 0 for j in 1:n for i in 1:m) &&
           all(b[i] >= 0 for i in 1:m) &&
           GBTools.has_slacks(A)
end

function parse_truncation(
    truncation_type::Symbol,
    A::Array{Int,2},
    b::Vector{Int}
)::Tuple{Symbol,DataType}
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

function groebner_basis(
    filepath :: String,
    markov_basis :: Union{Vector{Vector{Int}}, Nothing} = nothing;
    kwargs...
) :: Vector{Vector{Int}}
    optimizer = DEFAULT_SOLVER
    if :optimizer in keys(kwargs)
        optimizer = kwargs[:optimizer]
    end
    instance = IPInstance(filepath, optimizer=optimizer)
    return groebner_basis(instance, markov_basis; kwargs...)
end

function groebner_basis(
    model :: JuMP.Model,
    markov_basis :: Union{Vector{Vector{Int}}, Nothing} = nothing;
    kwargs...
) :: Vector{Vector{Int}}
    optimizer = DEFAULT_SOLVER
    if :optimizer in keys(kwargs)
        optimizer = kwargs[:optimizer]
    end
    instance = IPInstance(model, optimizer=optimizer)
    return groebner_basis(instance, markov_basis; kwargs...)
end

function groebner_basis(
    A :: Matrix{Int},
    b :: Vector{Int},
    C :: Matrix{T},
    markov_basis :: Union{Vector{Vector{Int}}, Nothing} = nothing;
    u :: Union{Nothing, Vector{<:Union{Int, Nothing}}} = nothing,
    kwargs...
) :: Vector{Vector{Int}} where {T <: Real}
    normalize = true
    if :apply_normalization in keys(kwargs)
        normalize = kwargs[:apply_normalization]
    end
    optimizer = DEFAULT_SOLVER
    if :optimizer in keys(kwargs)
        optimizer = kwargs[:optimizer]
    end
    instance = IPInstance(A, b, C, u, apply_normalization=normalize, optimizer=optimizer)
    return groebner_basis(instance, markov_basis; kwargs...)
end

@kwdef struct Options
    quiet :: Bool = true
    use_signatures :: Bool = false
    implicit_representation :: Bool = false
    truncation_type :: Symbol = :Heuristic
    use_quick_truncation :: Bool = true
    use_binary_truncation :: Bool = true
    solutions :: Vector{Vector{Int}} = Vector{Int}[]
    module_order :: Symbol = :ltpot
    auto_reduce_type :: Symbol = :FixedIterations
    auto_reduce_freq :: Float64 = 2500.0
    cache_tree_size :: Int = 100
    debug :: Bool = false
    info :: Bool = false
    time_limit :: Float64 = 0.0
    gb_size_limit :: Int = 0
    pair_processing :: Symbol = :FIFO
end

"""
    groebner_basis(
    instance :: IPInstance,
    markov_basis :: Union{Vector{Vector{Int}}, Nothing} = nothing;
    kwargs...
) :: Vector{Vector{Int}}

    Return a Gröbner basis of `instance` using the given `markov_basis` if
    it is available. Otherwise, a Markov basis is computed first.

    Keyword parameters:
    - quiet: Whether to print information about the computation
    - use_signatures: Whether to use signatures in the computation (work in progress)
    - implicit_representation: Whether to represent slack variables implicitly or not (as described in Thomas and Weismantel (1997))
    - truncation_type: The type of truncation to use. Can be :Heuristic, :LP, :IP, :Simple, or :None
    - use_quick_truncation: Whether to use grading-based truncation
    - use_binary_truncation: Whether to use the binary truncation criterion for binary variables
    - solutions: A vector with known feasible solutions to be optimized
    - module_order: The module order to use in a signature GB computation (work in progress)
"""
function groebner_basis(
    instance :: IPInstance,
    markov_basis :: Union{Vector{Vector{Int}}, Nothing} = nothing;
    kwargs...
) :: Vector{Vector{Int}}
    #Check whether the given Markov basis is in the kernel of instance.A
    if !isnothing(markov_basis) && !in_kernel(markov_basis, instance)
        throw(ArgumentError("markov_basis is not in the kernel of instance."))
    end
    opt = Options(; kwargs...)
    #If there is no Markov basis given, compute it
    if isnothing(markov_basis)
        markov_basis = Markov.markov_basis(instance, quiet=opt.quiet)
    end
    #Build the GBAlgorithm corresponding to the given parameters

    algorithm_type = opt.use_signatures ? SignatureAlgorithm : BuchbergerAlgorithm
    representation = opt.implicit_representation ? GradedBinomial : Binomial
    use_minimization = true
    trunc_type, trunc_var = parse_truncation(opt.truncation_type, instance.A, instance.b)

    #Signatures aren't currently supported with implicit representation
    if opt.use_signatures && opt.implicit_representation
        throw(ArgumentError("Signatures are not supported with implicit representation."))
    end

    #Run GB algorithm over the given instance
    if opt.use_signatures
        #TODO: Change the SignatureAlgorithm constructor, then send the whole
        #instance
        algorithm = algorithm_type(
            representation, instance.C, instance.A, instance.b, opt.module_order,
            truncate, use_minimization
        )
    else
        algorithm = algorithm_type(
            markov_basis, instance, opt.solutions, T = representation, truncation_type = trunc_type,
            trunc_var_type = trunc_var, minimization = use_minimization,
            use_quick_truncation = opt.use_quick_truncation, use_binary_truncation = opt.use_binary_truncation,
            auto_reduce_type = opt.auto_reduce_type, auto_reduce_freq = opt.auto_reduce_freq,
            cache_tree_size = opt.cache_tree_size, debug = opt.debug, info = opt.info,
            pair_processing = opt.pair_processing
        )
    end
    gb = GBAlgorithms.run(
        algorithm, quiet = opt.quiet, time_limit = opt.time_limit,
        gb_size_limit = opt.gb_size_limit
    )
    @debug "IPGBs finished, GB:" gb
    return gb
end

end
