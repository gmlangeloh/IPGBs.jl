module GBAlgorithms

export GBAlgorithm, run, current_basis, update!, stats, order, increment,
    truncate_basis

using JuMP

using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances

using IPGBs
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Orders

"""
A generic Gröbner Basis algorithm. Currently, it can be either a
BuchbergerAlgorithm or a SignatureAlgorithm.
"""
abstract type GBAlgorithm end

# GBAlgorithm interface
stats(algorithm :: GBAlgorithm) = algorithm.stats
elem(algorithm :: GBAlgorithm) = algorithm.preallocated #Memory allocated for current S-pair
next_pair!(:: GBAlgorithm) = error("Not implemented.")
current_basis(algorithm :: GBAlgorithm) = algorithm.basis
update!(:: GBAlgorithm, :: GBElement, :: Union{CriticalPair, Nothing}) =
    error("Not implemented.")
truncate_basis(algorithm :: GBAlgorithm) = algorithm.should_truncate
use_implicit_representation(:: GBAlgorithm) = false
quick_truncation(:: GBAlgorithm, :: GBElement) = false
optimize_solutions!(:: GBAlgorithm) = error("Not implemented.")

function initialize!(
    :: GBAlgorithm,
    :: Array{Int, 2},
    :: Vector{Int},
    :: Array{Int, 2},
    :: Vector{Int}
)
    error("Not implemented.")
end

"""
    increment(algorithm :: GBAlgorithm, stat :: Symbol)

Increment the field `stat` of stats(`algorithm`)
"""
function increment(
    algorithm :: GBAlgorithm,
    stat :: Symbol
)
    alg_stats = stats(algorithm)
    old_value = getfield(alg_stats, stat)
    setfield!(alg_stats, stat, old_value + 1)
end

#The following functions may or may not be extended, depending on the algorithm.
"""
    late_pair_elimination(:: GBAlgorithm, :: CriticalPair)

Return true iff the given CriticalPair should be eliminated according to this
algorithm's criteria.
"""
late_pair_elimination(:: GBAlgorithm, :: CriticalPair) = nothing

process_zero_reduction!(:: GBAlgorithm, :: T, :: CriticalPair) where {T <: GBElement} = nothing

function sbinomial(
    algorithm :: GBAlgorithm,
    pair :: CriticalPair
) :: GBElement
    increment(algorithm, :built_pairs)
    return BinomialSets.sbinomial(elem(algorithm), pair, current_basis(algorithm))
end

function reduce!(
    algorithm :: GBAlgorithm,
    binomial :: GBElement
) :: Tuple{Bool, Bool}
    increment(algorithm, :reduced_pairs)
    return BinomialSets.reduce!(binomial, current_basis(algorithm))
end

function truncate(
    algorithm :: GBAlgorithm,
    binomial :: GBElement
) :: Bool
    instance = algorithm.instance
    truncated = GBElements.truncate(
        binomial, instance.A, instance.b, instance.u, algorithm.model,
        algorithm.model_constraints, algorithm.should_truncate,
        algorithm.truncation_type
    )
    if truncated
        increment(algorithm, :eliminated_by_truncation)
        return true
    end
    return false
end

"""
    reintroduce_truncated!(algorithm :: GBAlgorithm)

Add truncated elements of the given Markov basis back into the GB.

This is used at the end of the computation to keep the GB as a basis of
the input ideal, regardless of truncation.
"""
function reintroduce_truncated!(
    algorithm :: GBAlgorithm
)
    for gen in algorithm.truncated_gens
        update!(algorithm, gen, nothing)
    end
end

"""
    is_minimization(algorithm :: GBAlgorithm) :: Bool

Return true iff this algorithm is working with problems in minimization form.
"""
function is_minimization(
    algorithm :: GBAlgorithm
) :: Bool
    return BinomialSets.is_minimization(current_basis(algorithm))
end

"""
    prepare_gb_output(algorithm :: GBAlgorithm) :: Vector{Vector{Int}}

Return the user-friendly output Gröbner Basis of `algorithm`, reduced,
sorted and in the original variable ordering.

This output is compatible with 4ti2's output format.
"""
function prepare_gb_output(
    algorithm :: GBAlgorithm
) :: Vector{Vector{Int}}
    reintroduce_truncated!(algorithm)
    reduced_basis!(current_basis(algorithm))
    optimize_solutions!(algorithm)
    output = fourti2_form(current_basis(algorithm))
    sort!(output)
    return output
end

function print_algorithm_stats(
    algorithm :: GBAlgorithm,
    quiet :: Bool
)
    @info "GB Algorithm stats:" stats(algorithm)
    if !quiet
        println(stats(algorithm))
        #This is kind of a hack but works
        #println(current_basis(algorithm).reduction_tree.stats)
    end
end

"""
    run(algorithm :: GBAlgorithm, quiet :: Bool) :: Vector{Vector{Int}}

Return the Gröbner Basis obtained by running `algorithm`.
"""
function run(
    algorithm :: GBAlgorithm;
    quiet :: Bool = false
) :: Vector{Vector{Int}}
    #Main loop: process all relevant S-pairs
    #TODO: Move parameter initialization elsewhere!!! (IPGBs.jl)
    IPGBs.initialize_parameters(
        auto_reduce_freq = 2500,
        auto_reduce_type = IPGBs.FIXED_ITERATIONS,
        cache_tree_size = 500,
        debug = false,
        info = false
    )
    while true
        pair = next_pair!(algorithm)
        if isnothing(pair) #All S-pairs were processed, terminate algorithm.
            break
        end
        if late_pair_elimination(algorithm, pair)
            continue
        end
        #Generate S-pair, reduce it and add to basis if necessary
        binomial = sbinomial(algorithm, pair)
        if quick_truncation(algorithm, binomial)
            continue
        end
        reduced_to_zero, _ = reduce!(algorithm, binomial)
        if !reduced_to_zero
            #First we check reduction to zero, then truncation.
            #This is more efficient because relatively few elements are
            #truncated when compared to the number of zero reductions
            if !truncate(algorithm, binomial)
                update!(algorithm, binomial, pair)
            end
        else
            process_zero_reduction!(algorithm, binomial, pair)
        end
    end
    print_algorithm_stats(algorithm, quiet)
    #tr(bin) = quick_truncation(algorithm, bin) || truncate(algorithm, bin)
    #ans = is_truncated_groebner_basis(current_basis(algorithm), tr)
    #println("Is this a truncated GB? " * string(ans))
    output = prepare_gb_output(algorithm)
    #println("nodes ", current_basis(algorithm).reduction_tree.stats.reduction_steps)
    return output
end

end
