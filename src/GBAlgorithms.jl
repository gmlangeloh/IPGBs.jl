module GBAlgorithms

export GBAlgorithm, run, current_basis, update!, GBStats, stats, order, increment,
    truncate_basis

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.Orders

abstract type GBStats end

function Base.show(
    io :: IO,
    stats :: GBStats
)
    println(io, "Algorithm statistics:")
    i = 1
    fields = fieldnames(typeof(stats))
    num_fields = length(fields)
    for field in fields
        if i < num_fields
            println(io, field, " => ", getfield(stats, field))
        else
            print(io, field, " => ", getfield(stats, field))
        end
        i += 1
    end
end

"""
A generic Gröbner Basis algorithm. For now, it can be either a BuchbergerAlgorithm
or a SignatureAlgorithm.
"""
abstract type GBAlgorithm end

# GBAlgorithm interface
stats(algorithm :: GBAlgorithm) = algorithm.stats
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
    return BinomialSets.sbinomial(pair, current_basis(algorithm))
end

function reduce!(
    algorithm :: GBAlgorithm,
    binomial :: GBElement
) :: Bool
    increment(algorithm, :reduced_pairs)
    return BinomialSets.reduce!(binomial, current_basis(algorithm))
end

function truncate(
    algorithm :: GBAlgorithm,
    binomial :: GBElement,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool
    if !truncate_basis(algorithm)
        return false
    end
    should_truncate = !isfeasible(binomial, A, b, u)
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
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    quiet :: Bool = false
) :: Vector{Vector{Int}}
    #Compute the initial basis of the toric ideal
    A, b, C, u = GBTools.normalize(
        A, b, C, u, apply_normalization=!use_implicit_representation(algorithm),
        invert_objective=is_minimization(algorithm)
    )
    initialize!(algorithm, A, b, C, u)
    #Main loop: process all relevant S-pairs
    while true
        pair = next_pair!(algorithm)
        if isnothing(pair) #All S-pairs were processed, terminate algorithm.
            break
        end
        #println(pair.j - 1, " ", pair.i - 1)
        #println("u = ", current_basis(algorithm)[pair.j])
        #println("v = ", current_basis(algorithm)[pair.i])
        if late_pair_elimination(algorithm, pair)
            #println("eliminated")
            continue
        end
        binomial = sbinomial(algorithm, pair)
        #if !truncate(algorithm, binomial, A, b, u)
        reduced_to_zero = reduce!(algorithm, binomial)
        if !reduced_to_zero && !truncate(algorithm, binomial, A, b, u)
            update!(algorithm, binomial, pair)
            #println("basis")
        elseif reduced_to_zero #Update syzygies in case this makes sense
            process_zero_reduction!(algorithm, binomial, pair)
            #println("zero")
        else
            #println("truncated")
        end
        #end
    end
    if !quiet
        println(stats(algorithm))
    end
    reduced_basis!(current_basis(algorithm))
    return fourti2_form(current_basis(algorithm))
end

end
