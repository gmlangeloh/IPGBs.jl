module BinomialSets

export BinomialSet, order, binomials, reduction_tree,
    is_support_reducible, fourti2_form, sbinomial, minimal_basis!, reduced_basis!,
    auto_reduce_once!, is_minimization, change_ordering!, is_groebner_basis,
    is_truncated_groebner_basis

using IPGBs.FastBitSets
using IPGBs.IPInstances
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

function BinomialSet(
    basis :: Vector{T},
    C :: Array{S},
    A :: Matrix{Int},
    b :: Vector{Int}
) where {T <: AbstractVector{Int}, S <: Real}
    is_minimization = !is_implicit(T)
    if !(S <: Float64)
        C = Float64.(C)
    end
    order = MonomialOrder(C, A, b, is_minimization)
    return BinomialSet{T, MonomialOrder}(basis, order, is_minimization)
end

"""
Changes the ordering of `bs` to a monomial order given by the matrix
`new_order`.
"""
function change_ordering!(
    bs :: BinomialSet{T, S},
    new_order :: Array{Float64, 2},
    A :: Array{Int, 2},
    b :: Vector{Int}
) where {T <: AbstractVector{Int}, S <: GBOrder}
    Orders.change_ordering!(bs.order, new_order, A, b)
    #TODO I probably should also reorientate all elements of bs
end

function binomials(fourti2_basis :: Matrix{Int})
    return [ fourti2_basis[i, :] for i in 1:size(fourti2_basis, 1) ]
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
Reduces `g` by `bs` efficiently using a Support Tree. Returns a pair of booleans
(r, c) where r == true iff `g` reduced to zero, and c == true iff `g` was
changed in this reduction.
"""
function reduce!(
    g :: T,
    bs :: BinomialSet{T, S};
    skipbinomial :: Union{T, Nothing} = nothing
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}, S <: GBOrder}
    return reduce!(g, bs, reduction_tree(bs), skipbinomial=skipbinomial)
end

"""
Fully reduce `binomial` by `gb` in-place, finding its normal form. Uses `tree`
to speed up the search for reducers. Returns a pair of booleans (r, c) where
r == true iff `g` reduced to zero, and c == true iff `g` was changed in this
reduction.

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
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    is_singular = Ref{Bool}(false)
    changed = false
    #Reduce both the leading and trailing terms, leading first
    for negative in (false, true)
        if negative && !is_minimization(gb)
            #In implicit form, the reduction already reduces both leading and
            #trailing term. In this case, there's nothing to do here.
            return false, changed
        end
        while true
            reducer = find_reducer(
                g, gb, tree, skipbinomial=skipbinomial, is_singular=is_singular,
                negative=negative
            )
            #g has a singular signature, so it reduces to zero
            #We can ignore the reducer and just say `g` reduces to zero
            #if haskey(params, "is_singular") && params["is_singular"]
            if is_singular[]
                return true, changed
            end
            #No reducer found, terminate search
            if isnothing(reducer)
                if negative
                    return false, changed
                end
                break
            end
            changed = true
            #Found some reducer, add it to histogram if one is available
            if !isnothing(reduction_count)
                for i in 1:length(gb)
                    if gb[i] === reducer
                        reduction_count[i] += 1
                        break
                    end
                end
            end
            if !GBElements.is_negative_disjoint(g, reducer, negative=negative)
                return true, changed
            end
            #Now apply the reduction and check if it is a zero reduction
            if has_order(gb)
                reduced_to_zero = GBElements.reduce!(
                    g, reducer, order(gb), negative=negative
                )
            else
                reduced_to_zero = GBElements.reduce!(
                    g, reducer, negative=negative
                )
            end
            if reduced_to_zero
                return true, changed
            end
        end
    end
    @assert false #We should never reach this
    return false, false
end

"""
Creates a concrete S-binomial from `pair`. In practice, this should only be
called after we were unable to eliminate `pair`.
"""
function sbinomial(
    mem :: Vector{Int},
    pair :: CriticalPair,
    bs :: BinomialSet{T, S}
) :: T where {T <: AbstractVector{Int}, S <: GBOrder}
    v = bs[GBElements.first(pair)]
    w = bs[GBElements.second(pair)]
    r = build(mem, v, w, pair)
    orientate!(r, order(bs))
    return r
end

"""
Returns true iff the i-th and j-th vectors of `bs` have disjoint negative
bounded components.

This is essentially 4ti2's cancellation criterion.
"""
function is_negative_disjoint(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S},
    instance :: IPInstance
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    u = bs[i]
    v = bs[j]
    for k in 1:instance.bounded_end
        if u[k] < 0 && v[k] < 0
            return false
        end
    end
    return true
end

"""
Returns true iff the i-th and j-th vectors of `bs` have disjoint positive
non-negative components.

This is essentially the GCD criterion.
"""
function is_positive_disjoint(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S},
    instance :: IPInstance
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    u = bs[i]
    v = bs[j]
    for k in 1:instance.nonnegative_end
        if u[k] > 0 && v[k] > 0
            return false
        end
    end
    return true
end

"""
Returns true if (i, j) should be discarded.

This is an implementation of Criteria 1 and 2 from Malkin's thesis.

Criterion 1 is simply the classical Buchberger GCD criterion. If the
positive supports are disjoint (= GCD of the leading terms is 1) then the
pair may be discarded.

Criterion 2 is specific to homogeneous ideals (or just ideals coming from
lattices?) Either way, if the negative supports are not disjoint (= GCD of
trailing terms is not 1) then the pair may be discarded. This applies
specifically to the bounded components (= variables) of the problem.
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S},
    instance :: IPInstance
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    neg_sup = bs.negative_supports
    pos_sup = bs.positive_supports
    #If there are unbounded components, do not use support bitsets
    if instance.bounded_end < instance.n
        if is_minimization(bs)
            return !is_negative_disjoint(i, j, bs, instance) ||
                is_positive_disjoint(i, j, bs, instance)
        end
        return is_negative_disjoint(i, j, bs, instance) ||
            !is_positive_disjoint(i, j, bs, instance)
    end
    #Use support bitsets for a little more efficiency
    if is_minimization(bs)
        return !disjoint(neg_sup[i], neg_sup[j]) ||
            disjoint(pos_sup[i], pos_sup[j])
    end
    #Maximization problem
    return disjoint(neg_sup[i], neg_sup[j]) ||
        !disjoint(pos_sup[i], pos_sup[j])
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
    mem = Vector{Int}(undef, length(bs[1]))
    for i in 1:length(bs)
        for j in 1:(i-1)
            s = BinomialPair(i, j)
            binomial = sbinomial(mem, s, bs)
            reduced_to_zero, _ = reduce!(binomial, bs)
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
    mem = Vector{Int}(undef, length(bs[1]))
    for i in 1:length(bs)
        for j in 1:(i-1)
            s = BinomialPair(i, j)
            binomial = sbinomial(mem, s, bs)
            reduced_to_zero, _ = reduce!(binomial, bs)
            #TODO we likely can't always use simple_truncation here.
            if !reduced_to_zero && simple_truncation(binomial, A, b, u)
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
"""
function auto_reduce_once!(
    gb :: BinomialSet{T, S},
    index :: Int
) :: Tuple{Int, Int} where {T <: AbstractVector{Int}, S <: GBOrder}
    removed = 0
    removed_before_index = 0
    for i in length(gb):-1:1
        g = gb[i]
        reduced_to_zero, changed = reduce!(g, gb, skipbinomial=g)
        if reduced_to_zero
            deleteat!(gb, i)
            removed += 1
            if i < index
                removed_before_index += 1
            end
        elseif changed
            #Change was already made in-place in reduce!
            #No need to update the basis in this case
        end
    end
    return removed, removed_before_index
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
    #Currently, this implementation doesn't support signatures
    #TODO Fix this later
    if has_signature(T)
        return
    end
    minimal_basis!(gb)
    if !is_minimization(gb)
        #Doesn't compute reduced basis for implicit representation
        #TODO fix this in some later version!
        return
    end
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
