module BinomialSets

export BinomialSet, order, binomials, reduction_tree,
    is_support_reducible, fourti2_form, sbinomial, minimal_basis!, reduced_basis!,
    auto_reduce!, is_minimization, change_ordering!, is_groebner_basis, is_truncated_groebner_basis

using IPGBs.FastBitSets
using IPGBs.Orders
using IPGBs.GBElements
using IPGBs.SupportTrees

"""
Represents a set of binomials (e.g. a partial or complete Gröbner Basis),
including structures for efficient reduction of elements.

In order to allow 4ti2-form bases, we allow elements to be Vector{Int},
in addition to GBElements.
"""
struct BinomialSet{T <: AbstractVector{Int}, S <: GBOrder} <: AbstractVector{T}
    basis :: Vector{T}
    order :: S
    reduction_tree :: SupportTree{T}
    minimization_form :: Bool #TODO maybe this should be part of the order?

    #We store the supports here instead of on the elements themselves to avoid
    #having to compute them unnecessarily or having to compute them after creating
    #elements and then updating these elements.
    positive_supports :: Vector{FastBitSet}
    negative_supports :: Vector{FastBitSet}

    function BinomialSet{T, S}(basis :: Vector{T}, order :: S, min :: Bool) where {T, S}
        tree = support_tree(basis, fullfilter=is_implicit(T))
        pos_supps, neg_supps = GBElements.supports(basis)
        new{T, S}(
            basis, order, tree, min, pos_supps, neg_supps
        )
    end
end

function BinomialSet(
    basis :: Vector{T},
    order :: S
) where {T <: AbstractVector{Int}, S <: GBOrder}
    return BinomialSet{T, S}(basis, order, !is_implicit(T))
end

"""
Creates a BinomialSet from `basis` with respect to `C`.

`module_order` is only relevant when `T` is a type of binomial with signature.
"""
function BinomialSet(
    basis :: Vector{T},
    C :: Array{Int, 2};
    module_order :: Symbol = :pot
) where {T <: AbstractVector{Int}}
    if has_signature(T)
        #TODO this doesn't work yet, due to cyclic dependencies between
        #this module and SignaturePolynomials. I need to find a way to
        #break some of these dependencies.
        #order = ModuleMonomialOrdering(C, module_order, basis)
    else
        order = MonomialOrder(C)
    end
    return BinomialSet{T, MonomialOrder}(basis, order)
end

"""
Changes the ordering of `bs` to a monomial order given by the matrix
`new_order`.
"""
function change_ordering!(
    bs :: BinomialSet{T, S},
    new_order :: Array{Int, 2}
) where {T <: AbstractVector{Int}, S <: GBOrder}
    Orders.change_ordering!(bs.order, new_order)
    #TODO I probably should also reorientate all elements of bs
end

has_order(:: BinomialSet) = true
has_order(:: AbstractVector) = false

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
) :: Tuple where {T <: AbstractVector{Int}, S <: GBOrder}
    return size(bs.basis)
end

function Base.getindex(
    bs :: BinomialSet{T, S},
    i :: Int
) :: T where {T <: AbstractVector{Int}, S <: GBOrder}
    return bs.basis[i]
end

function Base.length(
    bs :: BinomialSet{T, S}
) :: Int where {T <: AbstractVector{Int}, S <: GBOrder}
    return length(bs.basis)
end

function Base.push!(
    bs :: BinomialSet{T, S},
    g :: T
) where {T <: AbstractVector{Int}, S <: GBOrder}
    push!(bs.basis, g)
    addbinomial!(bs.reduction_tree, bs[length(bs)])
    p, n = GBElements.supports(g)
    push!(bs.positive_supports, p)
    push!(bs.negative_supports, n)
end

function Base.deleteat!(
    bs :: BinomialSet{T, S},
    i :: Int
) where {T <: AbstractVector{Int}, S <: GBOrder}
    deleted_elem = bs[i]
    removebinomial!(reduction_tree(bs), deleted_elem)
    deleteat!(bs.basis, i)
    deleteat!(bs.positive_supports, i)
    deleteat!(bs.negative_supports, i)
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
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    return reduce!(g, bs, reduction_tree(bs))
end

"""
Fully reduce `binomial` by `gb` in-place, finding its normal form. Uses `tree`
to speed up the search for reducers. Returns true iff `binomial` reduces to zero.

`binomial` can also be a monomial.

If `reduction_count` is given, the number of times each reducer was used
in this reduction process is added to reduction_count. This greatly slows
down execution, and is only meant for experimental purposes. The parameter
should be set to `nothing` in practical runs.
"""
function reduce!(
    g :: T,
    gb :: S,
    tree :: SupportTree{T};
    reduction_count :: Union{Vector{Int}, Nothing} = nothing,
    skipbinomial :: Union{T, Nothing} = nothing
) :: Bool where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    is_singular = Ref{Bool}(false)
    while true
        reducer = find_reducer(
            g, gb, tree, skipbinomial=skipbinomial, is_singular=is_singular
        )
        #g has a singular signature, so it reduces to zero
        #We can ignore the reducer and just say `g` reduces to zero
        #if haskey(params, "is_singular") && params["is_singular"]
        if is_singular[]
            return true
        end
        #No reducer found, terminate search
        if isnothing(reducer)
            return false
        end
        #Found some reducer, add it to histogram if one is available
        if !isnothing(reduction_count)
            for i in 1:length(gb)
                if gb[i] === reducer
                    reduction_count[i] += 1
                    break
                end
            end
        end
        #Now apply the reduction and check if it is a zero reduction
        if has_order(gb)
            reduced_to_zero = GBElements.reduce!(g, reducer, order(gb))
        else
            reduced_to_zero = GBElements.reduce!(g, reducer)
        end
        if reduced_to_zero
            return true
        end
    end
    @assert false #We should never reach this
    return false
end

"""
Creates a concrete S-binomial from `pair`. In practice, this should only be
called after we were unable to eliminate `pair`.
"""
function sbinomial(
    pair :: CriticalPair,
    bs :: BinomialSet{T, S}
) :: T where {T <: AbstractVector{Int}, S <: GBOrder}
    v = bs[GBElements.first(pair)]
    w = bs[GBElements.second(pair)]
    if Base.lt(order(bs), v, w)
        r = build(w, v, pair)
    else
        r = build(v, w, pair)
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
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
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
theorem and generating all S-binomials, thus may be slow. Note that this verifies if
`bs` is a Gröbner Basis, but not a truncated Gröbner Basis. To check the latter case,
use is_truncated_groebner_basis.

Useful for debugging and automatic tests.
"""
function is_groebner_basis(
    bs :: BinomialSet{T, S}
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
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
Checks whether `bs` is a truncated Gröbner Basis where truncation is done with
respect to

Ax <= b
x <= u
"""
function is_truncated_groebner_basis(
    bs :: BinomialSet{T, S},
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    for i in 1:length(bs)
        for j in 1:(i-1)
            s = BinomialPair(i, j)
            binomial = sbinomial(s, bs)
            reduced_to_zero = reduce!(binomial, bs)
            if !reduced_to_zero && isfeasible(binomial, A, b, u)
                #Note we check isfeasible after reduction. This is often
                #quicker, because then we skip the check for zero reductions
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
) :: Vector{Vector{Int}} where {T <: AbstractVector{Int}, S <: GBOrder}
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
) where {T <: AbstractVector{Int}, S <: GBOrder}
    for i in length(gb):-1:1
        nothing
    end
end

"""
Updates gb to a minimal Gröbner Basis.

For more info on minimal GBs, see Lemma 3 from Cox, Little, O'Shea Chapter 2.7.
In summary, it shows that one can remove from a GB any g such that LT(g) is a
multiple of LT(h), for h distinct from g in the GB.
"""
function minimal_basis!(
    gb :: BinomialSet{T, S}
) where {T <: AbstractVector{Int}, S <: GBOrder}
    for i in length(gb):-1:1
        g = gb[i]
        reducer = find_reducer(g, gb, reduction_tree(gb), skipbinomial=g)
        if !isnothing(reducer)
            deleteat!(gb, i)
        end
    end
end

"""
Updates gb to a reduced Gröbner Basis.
"""
function reduced_basis!(
    gb :: BinomialSet{T, S}
) where {T <: AbstractVector{Int}, S <: GBOrder}
    minimal_basis!(gb)
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        while reducing
            h = find_reducer(g, gb, reduction_tree(gb), negative=true)
            if !isnothing(h)
                GBElements.reduce!(g, h, order(gb), negative=true)
            else
                reducing = false
            end
        end
    end
end

end
