module BinomialSets

export GBOrder, BinomialSet, order, binomials, reduction_tree, is_support_reducible

using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.GradedBinomials
using IPGBs.SupportTrees

"""
This is specialized by MonomialOrder in case of Buchberger's algorithm and by
ModuleMonomialOrder in case of Signature-based algorithms.
"""
abstract type GBOrder <: Base.Ordering end

struct MonomialOrder <: GBOrder
    cost_matrix :: Array{Int, 2}
    #Probably should carry a tiebreaker around too
end

#TODO this can be in a separated file for utility functions or something
function supports(
    B :: Vector{T}
) :: Tuple{Vector{FastBitSet}, Vector{FastBitSet}} where {T <: GBElement}
    pos_supps = FastBitSet[]
    neg_supps = FastBitSet[]
    for g in B
        p, n = GBElements.supports(g)
        push!(pos_supps, p)
        push!(neg_supps, n)
    end
    return pos_supps, neg_supps
end

struct BinomialSet{T <: GBElement, S <: GBOrder} <: AbstractVector{T}
    basis :: Vector{T}
    order :: S
    reduction_tree :: SupportTree{T}

    #We store the supports here instead of on the elements themselves to avoid
    #having to compute them unnecessarily or having to compute them after creating
    #elements and then updating these elements.
    positive_supports :: Vector{FastBitSet}
    negative_supports :: Vector{FastBitSet}
    function BinomialSet(basis :: Vector{T}, order :: S) where {T, S}
        #TODO also try to build the order
        tree = support_tree(basis, fullfilter=(T == GradedBinomial))
        pos_supps, neg_supps = supports(basis)
        new{T, S}(
            basis, order, tree, pos_supps, neg_supps
        )
    end
end

#
# Accessors
#

binomials(bs :: BinomialSet) = bs.basis
order(bs :: BinomialSet) = bs.order
reduction_tree(bs :: BinomialSet) = bs.reduction_tree
#TODO Do I even need to access the others externally?

#
# Implementing the AbstractVector interface
#

function Base.size(
    bs :: BinomialSet{T, S}
) :: Tuple where {T <: GBElement, S <: GBOrder}
    return size(bs.basis)
end

function Base.getindex(
    bs :: BinomialSet{T, S},
    i :: Int
) :: T where {T <: GBElement, S <: GBOrder}
    return bs.basis[i]
end

function Base.length(
    bs :: BinomialSet{T, S}
) :: Int where {T <: GBElement, S <: GBOrder}
    return length(bs.basis)
end

function Base.push!(
    bs :: BinomialSet{T, S},
    g :: T
) where {T <: GBElement, S <: GBOrder}
    push!(bs.basis, g)
    #TODO update this! Not really compatible with a SigBasis yet
    #The monomial orders should have a trait has_generators or something
    #Although the array in there should be just a reference to the one I update
    #I tested this in the REPL, and I think it works
    #push!(gb.order.generators, g)
    addbinomial!(bs.reduction_tree, bs[length(bs)])
    p, n = GBElements.supports(g)
    push!(bs.positive_supports, p)
    push!(bs.negative_supports, n)
end

#
# Additional functions defined on a BinomialSet
#

"""
Returns true if (i, j) should be discarded.

In a maximization problem, if (i, j) ..

TODO actually document this thing, it's not that trivial
BTW, this name is terrible. It should make clear that this is the gcd criterion
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S},
    minimization :: Bool #Maybe this one should be part of the BinomialSet struct
) :: Bool where {T <: GBElement, S <: GBOrder}
    if minimization
        return !disjoint(bs.negative_supports[i], bs.negative_supports[j]) ||
            disjoint(bs.positive_supports[i], bs.positive_supports[j])
    end
    #Maximization problem
    return disjoint(bs.negative_supports[i], bs.negative_supports[j]) ||
        !disjoint(bs.positive_supports[i], bs.positive_supports[j])
end

"""
Returns true iff `bs` is a Gröbner Basis. This is checked by applying Buchberger's
theorem and generating all S-binomials, thus may be slow.

Useful for debugging and automatic tests.
"""
function is_groebner_basis(
    bs :: BinomialSet
) :: Bool
    #TODO implement something to generate S-pairs and do the Buchberger check

    return false
end

end
