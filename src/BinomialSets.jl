module BinomialSets

export GBOrder, BinomialSet, MonomialOrder, order, binomials, reduction_tree,
    is_support_reducible, fourti2_form, sbinomial, minimal_basis!, reduced_basis!,
    auto_reduce!, is_minimization

using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.SupportTrees

"""
This is specialized by MonomialOrder in case of Buchberger's algorithm and by
ModuleMonomialOrder in case of Signature-based algorithms.
"""
abstract type GBOrder <: Base.Ordering end

#TODO for now this is useless. It should probably do SOMETHING
#Like tiebreaking and helping orientate elements
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
    minimization_form :: Bool #TODO maybe this should be part of the order?

    #We store the supports here instead of on the elements themselves to avoid
    #having to compute them unnecessarily or having to compute them after creating
    #elements and then updating these elements.
    positive_supports :: Vector{FastBitSet}
    negative_supports :: Vector{FastBitSet}

    function BinomialSet(basis :: Vector{T}, order :: S) where {T, S}
        #TODO also try to build the order
        tree = support_tree(basis, fullfilter=is_implicit(T))
        pos_supps, neg_supps = supports(basis)
        min = !is_implicit(T)
        new{T, S}(
            basis, order, tree, min, pos_supps, neg_supps
        )
    end
end

function BinomialSet(T :: Type)
    #TODO define order and minimization somewhere else, maybe pass them to this
    #constructor
    return BinomialSet{T, S}(T[])
end

#
# Accessors
#

binomials(bs :: BinomialSet) = bs.basis
order(bs :: BinomialSet) = bs.order
reduction_tree(bs :: BinomialSet) = bs.reduction_tree
is_minimization(bs :: BinomialSet) = bs.minimization_form

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
    addbinomial!(bs.reduction_tree, bs[length(bs)])
    p, n = GBElements.supports(g)
    push!(bs.positive_supports, p)
    push!(bs.negative_supports, n)
end

#
# Additional functions defined on a BinomialSet
#

"""
Reduces `g` by `bs` efficiently using a Support Tree. Returns true iff `g`
reduced to zero.
"""
function reduce!(
    g :: T,
    bs :: BinomialSet{T, S}
) :: Bool where {T <: GBElement, S <: GBOrder}
    return SupportTrees.reduce!(g, bs, reduction_tree(bs))
end

"""
Creates a concrete S-binomial from `pair`. In practice, this should only be
called after we were unable to eliminate `pair`.
"""
function sbinomial(
    pair :: CriticalPair,
    bs :: BinomialSet{T, S}
) :: T where {T <: GBElement, S <: GBOrder}
    v = bs[GBElements.first(pair)]
    w = bs[GBElements.second(pair)]
    #TODO maybe this complicated logic should be written somewhere else in terms
    #of orientate! or something similar?
    if cost(v) < cost(w)
        r = build(w, v, pair)
    elseif cost(w) < cost(v)
        r = build(v, w, pair)
    else #w.cost == v.cost
        if GBElements.lt_tiebreaker(v, w)
            r = build(w, v, pair)
        else
            r = build(v, w, pair)
        end
    end
    return r
end

"""
Returns true if (i, j) should be discarded.

In a maximization problem, if (i, j) ..

TODO actually document this thing, it's not that trivial
BTW, this name is terrible. It should make clear that this is the gcd criterion
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S}
) :: Bool where {T <: GBElement, S <: GBOrder}
    if is_minimization(bs)
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
    bs :: BinomialSet{T, S}
) :: Bool where {T <: GBElement, S <: GBOrder}
    for i in 1:length(bs)
        for j in 1:(i-1)
            s = BinomialPair(i, j)
            binomial = sbinomial(s, bs)
            reduced_to_zero = reduce!(binomial, bs)
            if !reduced_to_zero
                return false
            end
        end
    end
    return true
end

"""
Returns the binomials of `bs` as integer vectors in minimization form.
"""
function fourti2_form(
    bs :: BinomialSet{T, S}
) :: Vector{Vector{Int}} where {T <: GBElement, S <: GBOrder}
    sign = is_minimization(bs) ? 1 : -1
    return [ sign * fullform(g) for g in bs ]
end

"""
Reduces each element of the GB by the previous elements.

TODO count number of removed elements so we can decrement the iteration
index in the main loop
"""
function auto_reduce!(
    gb :: BinomialSet{T, S}
) where {T <: GBElement, S <: GBOrder}
    for i in length(gb):-1:1
        nothing
    end
end

"""
Updates gb to a minimal Gröbner Basis.

For more info on minimal GBs, see Lemma 3 from Cox, Little, O'Shea Chapter 2.7.
In summary, it shows that one can remove from a GB any g such that LT(g) is a
multiple of LT(h), for h distinct from g in the GB.

TODO doesn't work because deleteat! is not defined.
"""
function minimal_basis!(
    gb :: BinomialSet{T, S}
) where {T <: GBElement, S <: GBOrder}
    for i in length(gb):-1:1
        g = gb[i]
        reducer = find_reducer(g, gb, reduction_tree(gb), skipbinomial=g)
        if !isnothing(reducer)
            @show reducer
            @show g
            deleteat!(gb, i)
            removebinomial!(reduction_tree(gb), g)
        end
    end
end

"""
Updates gb to a reduced Gröbner Basis.

TODO bugfix this, it doesn't terminate or is just buggy overall
"""
function reduced_basis!(
    gb :: BinomialSet{T, S}
) where {T <: GBElement, S <: GBOrder}
    minimal_basis!(gb)
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        while reducing
            h = find_reducer(g, gb, reduction_tree(gb), negative=true)
            if !isnothing(h)
                GBElements.reduce!(g, h, negative=true)
            else
                reducing = false
            end
        end
    end
end

end
