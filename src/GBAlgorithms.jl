module GBAlgorithms

export GBAlgorithm, run, current_basis, update!, stats, order, increment,
    truncate_basis

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.IPInstances
using IPGBs.Orders

"""
A generic Gröbner Basis algorithm. For now, it can be either a BuchbergerAlgorithm
or a SignatureAlgorithm.
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
Increments the field `stat` of stats(algorithm)
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
Returns true iff the given CriticalPair should be eliminated according to this
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
    A = algorithm.instance.A
    b = algorithm.instance.b
    u = algorithm.instance.u
    if !truncate_basis(algorithm)
        return false
    end
    if algorithm.truncation_type == :Simple
        should_truncate = !simple_truncation(binomial, A, b, u)
    else #algorithm.truncation_type == :Model
        should_truncate = !model_truncation(
            binomial, A, b, algorithm.model, algorithm.model_constraints
        )
    end
    if should_truncate
        increment(algorithm, :eliminated_by_truncation)
        return true
    end
    return false
end

"""
Returns true iff this algorithm is working with problems in minimization form.
"""
function is_minimization(
    algorithm :: GBAlgorithm
) :: Bool
    return BinomialSets.is_minimization(current_basis(algorithm))
end

# Main GB algorithm logic. Does not depend on the specific algorithm.
function run(
    algorithm :: GBAlgorithm,
    instance :: IPInstance;
    quiet :: Bool = false
) :: Vector{Vector{Int}}
    #TODO move the computation of the initial basis somewhere else
    initialize!(algorithm)
    #Main loop: process all relevant S-pairs
    while true
        pair = next_pair!(algorithm)
        if isnothing(pair) #All S-pairs were processed, terminate algorithm.
            break
        end
        if late_pair_elimination(algorithm, pair)
            continue
        end
        binomial = sbinomial(algorithm, pair)
        reduced_to_zero, _ = reduce!(algorithm, binomial)
        if !reduced_to_zero && !truncate(algorithm, binomial)
            update!(algorithm, binomial, pair)
        elseif reduced_to_zero #Update syzygies in case this makes sense
            process_zero_reduction!(algorithm, binomial, pair)
        end
    end
    if !quiet
        println(stats(algorithm))
        #This is kind of a hack but works
        println(current_basis(algorithm).reduction_tree.stats)
    end
    reduced_basis!(current_basis(algorithm))
    output = fourti2_form(current_basis(algorithm))
    sort!(output)
    #TODO We need to reintroduce generators which were truncated!
    #To do this, I need to change some stuff in initialize / lattice_generator
    return output
end

end
