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
next_sbinomial(:: GBAlgorithm) = error("Not implemented.")
current_basis(:: GBAlgorithm) = error("Not implemented.")

#These may or may not be extended, depending on the algorithm.
process_zero_reduction(:: GBAlgorithm, :: T) where {T <: GBElement} = nothing

# Main GB algorithm logic. Does not depend on the specific algorithm.
function run(
    algorithm :: GBAlgorithm,
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    gb = current_basis(algorithm)
    while true
        #TODO need to consider all kinds of elimination criteria here
        s = next_sbinomial(algorithm)
        if isnothing(s)
            break
        end
        if !isfeasible(s, A, b, u)
            continue
        end
        reduced_to_zero = BinomialSets.reduce!(s, gb)
        if !reduced_to_zero
            push!(gb, s)
        else #Update syzygies in case this makes sense
            process_zero_reduction(algorithm)
        end
    end
    minimal_basis!(gb)
    return fourti2_form(gb)
end

end
