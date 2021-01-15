module GBAlgorithms

export GBAlgorithm, CriticalPair, run, current_basis

using IPGBs.BinomialSets
using IPGBs.GBElements

abstract type CriticalPair end

#CriticalPair interface
first(:: CriticalPair) :: Int = error("Not implemented.")
second(:: CriticalPair) :: Int = error("Not implemented.")

"""
Creates a concrete S-binomial from `pair`. In practice, this should only be
called after we were unable to eliminate `pair`.
"""
function sbinomial(
    pair :: CriticalPair,
    bs :: BinomialSet{T, S}
) :: T where {T <: GBElement, S <: GBOrder}
    #Probably won't work for Signatures. Lacks the signature in the construction
    v = bs[first(pair)]
    w = bs[second(pair)]
    if cost(v) < cost(w)
        r = w - v #TODO these are relatively expensive
    elseif cost(w) < cost(v)
        r = v - w
    else #w.cost == v.cost
        if GBElements.lt_tiebreaker(v, w)
            r = w - v
        else
            r = v - w
        end
    end
    return r
end

"""
A generic Gröbner Basis algorithm. For now, it can be either a BuchbergerAlgorithm
or a SignatureAlgorithm.
"""
abstract type GBAlgorithm end

# GBAlgorithm interface
next_pair(:: GBAlgorithm) = error("Not implemented.")
current_basis(:: GBAlgorithm) = error("Not implemented.")

#The following functions may or may not be extended, depending on the algorithm.
"""
Returns true iff the given CriticalPair should be eliminated according to this
algorithm's criteria.
"""
late_pair_elimination(:: GBAlgorithm, :: CriticalPair) = nothing

process_zero_reduction(:: GBAlgorithm, :: T) where {T <: GBElement} = nothing

# Main GB algorithm logic. Does not depend on the specific algorithm.
function run(
    algorithm :: GBAlgorithm,
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    #TODO need complicated initialization around here!
    gb = current_basis(algorithm)
    while true
        pair = next_pair(algorithm)
        if isnothing(pair) #All SPairs were processed, terminate algorithm.
            break
        end
        #TODO cannot forget to deal with repeated signatures here!
        if late_pair_elimination(algorithm, pair)
            continue
        end
        binomial = sbinomial(pair, gb)
        if isfeasible(binomial, A, b, u)
            reduced_to_zero = BinomialSets.reduce!(binomial, gb)
            if !reduced_to_zero
                push!(gb, binomial)
                #TODO update pairs and other structures!!!
            else #Update syzygies in case this makes sense
                process_zero_reduction(algorithm, binomial)
            end
        end
    end
    minimal_basis!(gb)
    return fourti2_form(gb)
end

end
