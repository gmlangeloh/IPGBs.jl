#TODO: re-export the relevant functions: the 4ti2 interface and my GB
#implementations
module IPGBs
export groebner_basis

include("./Globals.jl")
include("./FastBitSets.jl")
include("./ZDDs.jl")
include("./FastestBitSets.jl")
include("./cbitsets/CBitSets.jl")
include("./Statistics.jl")
include("./Orders.jl")
include("./MonomialHeaps.jl")
include("./FastComparator.jl")
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

include("./FGLM.jl")
include("./Buchberger.jl")
include("./SignatureAlgorithms.jl")
include("./Markov.jl")
include("./StandardDecomposition.jl")

using .GBAlgorithms
using .Markov
import .Buchberger: BuchbergerAlgorithm
import .SignatureAlgorithms: SignatureAlgorithm
import .Binomials: Binomial
import .GradedBinomials: GradedBinomial

using JuMP
using Random
using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances

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

"""
Computes a Markov basis for `instance` and then uses it to compute a Gröbner
basis.
"""
function groebner_basis(
    A::Array{Int,2},
    b::Vector{Int},
    C::Array{T,2};
    u::Union{Nothing,Vector{<:Union{Int,Nothing}}} = nothing,
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true,
    normalize_ip::Bool = true
)::Vector{Vector{Int}} where {T<:Real}
    if isnothing(u)
        u = Union{Int,Nothing}[]
        for _ in 1:size(A, 2)
            push!(u, nothing)
        end
    end
    initialize_parameters(
        auto_reduce_freq = 2500,
        debug = false,
        info = false
    )
    instance = IPInstance(A, b, C, u, apply_normalization=normalize_ip)
    markov = markov_basis(instance)
    return groebner_basis(
        markov, instance,
        use_signatures = use_signatures,
        implicit_representation = implicit_representation,
        module_order = module_order,
        truncation_type = truncation_type,
        quiet = quiet
    )
end

"""
Computes a Markov basis for `instance` and then uses it to compute a Gröbner
basis.
"""
function groebner_basis(
    instance::IPInstance;
    solutions :: Vector{Vector{Int}} = Vector{Int}[],
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}}
    if !isempty(solutions)
        markov = markov_basis(instance, solution=solutions[1])
    else
        markov = markov_basis(instance)
    end
    return groebner_basis(
        markov, instance, solutions,
        use_signatures = use_signatures,
        implicit_representation = implicit_representation,
        module_order = module_order,
        truncation_type = truncation_type,
        quiet = quiet
    )
end

"""
Computes a Markov basis for `instance` and then uses it to compute a Gröbner
basis.
"""
function groebner_basis(
    model::JuMP.Model;
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}}
    return groebner_basis(
        IPInstance(model),
        use_signatures = use_signatures,
        implicit_representation = implicit_representation,
        module_order = module_order,
        truncation_type = truncation_type,
        quiet = quiet
    )
end

function groebner_basis(
    path :: String;
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}}
    return groebner_basis(
        IPInstance(path),
        use_signatures = use_signatures,
        implicit_representation = implicit_representation,
        module_order = module_order,
        truncation_type = truncation_type,
        quiet = quiet
    )
end

function groebner_basis(
    markov_basis::Vector{Vector{Int}},
    A::Array{Int,2},
    b::Vector{Int},
    C::Array{T,2};
    u::Union{Nothing,Vector{<:Union{Int,Nothing}}} = nothing,
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}} where {T<:Real}
    if isnothing(u)
        u = Union{Int,Nothing}[]
        for _ in 1:size(A, 2)
            push!(u, nothing)
        end
    end
    instance = IPInstance(A, b, C, u)
    return groebner_basis(
        markov_basis, instance,
        use_signatures = use_signatures,
        implicit_representation = implicit_representation,
        module_order = module_order,
        truncation_type = truncation_type,
        quiet = quiet
    )
end

"""
    groebner_basis(
    markov_basis::Vector{Vector{Int}},
    instance::IPInstance,
    solutions::Vector{Vector{Int}} = Vector{Int}[];
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}}

Return the reduced Gröbner Basis of the ideal generated by `markov`,
corresponding to the given integer programming instance. If `solutions`
are given, they are reduced by the Gröbner basis before being returned.

- use_signatures: Whether to use signatures for S-binomial elimination. This
is currently not finished.
- implicit_representation: Whether to use implicit variables for slacks and
upper bounds. This is currently unsupported.
- module_order: The options are :ltpot, :pot and :top. Only relevant for signature
algorithms.
- truncation_type: The options are :None, :Heuristic, :LP, :IP, :Simple.
- quiet: Whether to print debug information.
"""
function groebner_basis(
    markov_basis::Vector{Vector{Int}},
    instance::IPInstance,
    solutions::Vector{Vector{Int}} = Vector{Int}[];
    use_signatures::Bool = false,
    implicit_representation::Bool = false,
    module_order::Symbol = :ltpot,
    truncation_type::Symbol = :Heuristic,
    quiet::Bool = true
)::Vector{Vector{Int}}
    #Setting parameters
    @debug "Starting to compute Gröbner basis for: " markov_basis
    algorithm_type = use_signatures ? SignatureAlgorithm : BuchbergerAlgorithm
    representation = implicit_representation ? GradedBinomial : Binomial
    use_minimization = true
    trunc_type, trunc_var = parse_truncation(truncation_type, instance.A, instance.b)

    #Signatures aren't currently supported with implicit representation
    @assert !(use_signatures && implicit_representation)

    #Run GB algorithm over the given instance
    if use_signatures
        #TODO: Change the SignatureAlgorithm constructor, then send the whole
        #instance
        algorithm = algorithm_type(
            representation, instance.C, instance.A, instance.b, module_order,
            truncate, use_minimization
        )
    else
        algorithm = algorithm_type(
            markov_basis, instance, solutions, T = representation, truncation_type = trunc_type,
            trunc_var_type = trunc_var, minimization = use_minimization
        )
    end
    gb = GBAlgorithms.run(algorithm, quiet = quiet)
    @debug "IPGBs finished, GB:" gb
    return gb
end

end
