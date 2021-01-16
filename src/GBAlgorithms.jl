module GBAlgorithms

export GBAlgorithm, run, current_basis

using IPGBs.BinomialSets
using IPGBs.GBElements

"""
A generic Gröbner Basis algorithm. For now, it can be either a BuchbergerAlgorithm
or a SignatureAlgorithm.
"""
abstract type GBAlgorithm end

# GBAlgorithm interface
next_pair(:: GBAlgorithm) = error("Not implemented.")
current_basis(:: GBAlgorithm) = error("Not implemented.")
initialize_order(:: GBAlgorithm, :: Array{Int, 2}) = error("Not implemented")

#The following functions may or may not be extended, depending on the algorithm.
"""
Returns true iff the given CriticalPair should be eliminated according to this
algorithm's criteria.
"""
late_pair_elimination(:: GBAlgorithm, :: CriticalPair) = nothing

process_zero_reduction(:: GBAlgorithm, :: T) where {T <: GBElement} = nothing

function update!(
    algorithm :: GBAlgorithm,
    g :: T
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
end

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
                update!(algorithm, binomial)
            else #Update syzygies in case this makes sense
                process_zero_reduction(algorithm, binomial)
            end
        end
    end
    minimal_basis!(gb) #TODO I don't think this works with Signatures yet
    return fourti2_form(gb)
end

end
