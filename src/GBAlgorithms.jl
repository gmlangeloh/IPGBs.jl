module GBAlgorithms

export GBAlgorithm, run, current_basis, update!, GBStats, stats, data

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.GBTools

struct GBStats
    stats :: Dict{String, Any}
    function GBStats()
        d = Dict{String, Any}()
        #Initialize some stats used by every algorithm
        d["zero_reductions"] = 0
        d["max_basis_size"] = 0
        d["queued_pairs"] = 0
        d["built_pairs"] = 0
        d["reduced_pairs"] = 0
        d["eliminated_by_truncation"] = 0
        new(d)
    end
end

function Base.show(
    io :: IO,
    stats :: GBStats
)
    println(io, "Algorithm statistics:")
    i = 1
    for key in keys(stats.stats)
        if i < length(stats.stats)
            println(io, key, " => ", stats.stats[key])
        else
            print(io, key, " => ", stats.stats[key])
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
stats(:: GBAlgorithm) = error("Not implemented.") #Returns the GBStats object
data(:: GBAlgorithm) = error("Not implemented.") #Returns the dict GBStats wraps
next_pair!(:: GBAlgorithm) = error("Not implemented.")
current_basis(:: GBAlgorithm) = error("Not implemented.")
update!(:: GBAlgorithm, :: GBElement) = error("Not implemented.")

function initialize!(
    :: GBAlgorithm,
    :: Array{Int, 2},
    :: Vector{Int},
    :: Array{Int, 2},
    :: Vector{Int}
)
    error("Not implemented.")
end

#The following functions may or may not be extended, depending on the algorithm.
"""
Returns true iff the given CriticalPair should be eliminated according to this
algorithm's criteria.
"""
late_pair_elimination(:: GBAlgorithm, :: CriticalPair) = nothing

process_zero_reduction!(:: GBAlgorithm, :: T) where {T <: GBElement} = nothing

function sbinomial(
    algorithm :: GBAlgorithm,
    pair :: CriticalPair
) :: GBElement
    data(algorithm)["built_pairs"] += 1
    return BinomialSets.sbinomial(pair, current_basis(algorithm))
end

function reduce!(
    algorithm :: GBAlgorithm,
    binomial :: GBElement
) :: Bool
    data(algorithm)["reduced_pairs"] += 1
    return BinomialSets.reduce!(binomial, current_basis(algorithm))
end

function truncate(
    algorithm :: GBAlgorithm,
    binomial :: GBElement,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool
    should_truncate = !isfeasible(binomial, A, b, u)
    if should_truncate
        data(algorithm)["eliminated_by_truncation"] += 1
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
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    #Compute the initial basis of the toric ideal
    A, b, C, u = GBTools.normalize(
        A, b, C, u, apply_normalization=is_minimization(algorithm)
    )
    initialize!(algorithm, A, b, C, u)
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
        if !truncate(algorithm, binomial, A, b, u)
            reduced_to_zero = reduce!(algorithm, binomial)
            if !reduced_to_zero
                update!(algorithm, binomial)
            else #Update syzygies in case this makes sense
                process_zero_reduction!(algorithm, binomial)
            end
        end
    end
    println(stats(algorithm))
    #minimal_basis!(current_basis(algorithm)) #TODO I don't think this works with Signatures yet
    return fourti2_form(current_basis(algorithm))
end

end
