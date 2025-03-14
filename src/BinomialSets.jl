module BinomialSets

export BinomialSet, order, binomials, reduction_tree,
    is_support_reducible, fourti2_form, sbinomial, minimal_basis!, reduced_basis!,
    auto_reduce_once!, is_minimization, change_ordering!, is_groebner_basis,
    is_truncated_groebner_basis

using IPGBs.IPInstances
using IPGBs.Orders
using IPGBs.GBElements
using IPGBs.SupportTrees

using IPGBs.Binomials

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
    minimization_form :: Bool #TODO: maybe this should be part of the order?

    function BinomialSet{T, S}(basis :: Vector{T}, order :: S, min :: Bool) where {T, S}
        tree = SupportTree{T}(basis, fullfilter=is_implicit(T))
        GBElements.compute_supports.(basis)
        GBElements.compute_binaries.(basis)
        new{T, S}(basis, order, tree, min)
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
    ip :: IPInstance,
    S :: DataType = T
) where {T <: AbstractVector{Int}}
    is_minimization = !is_implicit(T)
    if isnothing(unbounded_variables(ip))
        unbounded = fill(false, ip.n)
    else
        unbounded = unbounded_variables(ip)
    end
    order = MonomialOrder(
        ip.C, ip.A, ip.b, unbounded, is_minimization, optimizer=ip.optimizer
    )
    oriented_basis = T[]
    initialize_binomials(ip, order)
    for g in basis
        oriented_g = copy(g)
        orientate!(oriented_g, order)
        push!(oriented_basis, oriented_g)
    end
    if S == Binomial
        new_oriented = [to_gbelement(v, order, S, false) for v in oriented_basis]
        return BinomialSet{S, MonomialOrder}(new_oriented, order, is_minimization)
    end
    return BinomialSet{T, MonomialOrder}(oriented_basis, order, is_minimization)
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
    #TODO: I probably should also reorientate all elements of bs
end

function binomials(fourti2_basis :: Matrix{Int})
    return [ fourti2_basis[i, :] for i in 1:size(fourti2_basis, 1) ]
end

function inverse_set(
    bs :: BinomialSet{T, S}
) :: BinomialSet{T, S} where {T <: AbstractVector{Int}, S <: GBOrder}
    inverted_elements = [ copy(element(g)) for g in bs.basis ]
    inverted_order = inverse_order(bs.order)
    inverted_basis = [ to_gbelement(v, inverted_order, T) for v in inverted_elements ]
    return BinomialSet{T, S}(inverted_basis, inverted_order, bs.minimization_form)
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
    add_binomial!(bs.reduction_tree, bs[length(bs)])
    GBElements.compute_supports(g)
    GBElements.compute_binaries(g)
end

function Base.insert!(
    bs :: BinomialSet{T, S},
    i :: Int,
    g :: T
) where {T <: AbstractVector{Int}, S <: GBOrder}
    insert!(bs.basis, i, g)
    add_binomial!(bs.reduction_tree, g)
    GBElements.compute_supports(g)
    GBElements.compute_binaries(g)
end

function Base.deleteat!(
    bs :: BinomialSet{T, S},
    i :: Int
) where {T <: AbstractVector{Int}, S <: GBOrder}
    deleted_elem = bs[i]
    remove_binomial!(reduction_tree(bs), deleted_elem)
    deleteat!(bs.basis, i)
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
    skipbinomial :: Union{T, Nothing} = nothing,
    is_monomial_reduction :: Bool = is_monomial(g)
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}, S <: GBOrder}
    return reduce!(
        g, bs, reduction_tree(bs), skipbinomial=skipbinomial,
        is_monomial_reduction=is_monomial_reduction
    )
end

"""
    reduce!(
    g :: T,
    bs :: Vector{T}
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}}

    Reduce `g` by `bs` inefficiently by linearly searching for a reducer.
    This is useful for very small examples and debugging. For anything else,
    consider making a BinomialSet and reducing using a SupportTree.

    Return a tuple (r, c) where r == true iff `g` reduced to zero, and c == true
    iff `g` was changed in this reduction.
"""
function reduce!(
    g :: T,
    bs :: Vector{T}
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}}
    function linear_reducer_search(g, bs)
        found_reducer = false
        reducer = copy(g)
        for r in bs
            if GBElements.reduces(r, g)
                found_reducer = true
                reducer = r
                break
            end
        end
        return reducer, found_reducer
    end
    reducer, found_reducer = linear_reducer_search(g, bs)
    changed = false
    reduced_to_zero = false
    while found_reducer
        changed = true
        GBElements.reduce!(g, reducer)
        if is_zero(g)
            reduced_to_zero = true
            break
        end
        reducer, found_reducer = linear_reducer_search(g, bs)
    end
    return reduced_to_zero, changed
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
    tree :: ReductionTree{T};
    reduction_count :: Union{Vector{Int}, Nothing} = nothing,
    skipbinomial :: Union{T, Nothing} = nothing,
    is_monomial_reduction :: Bool = false
) :: Tuple{Bool, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    changed = false
    #Reduce both the leading and trailing terms, leading first
    for negative in (false, true)
        if negative && !is_minimization(gb)
            #In implicit form, the reduction already reduces both leading and
            #trailing term. In this case, there's nothing to do here.
            return false, changed
        end
        while true
            reducer, found_reducer = find_reducer(
                g, gb, tree, skipbinomial=skipbinomial,
                negative=negative
            )
            #g has a singular signature, so it reduces to zero
            #We can ignore the reducer and just say `g` reduces to zero
            #if haskey(params, "is_singular") && params["is_singular"]
            #if is_singular[]
            #    return true, changed
            #end
            #No reducer found, terminate search
            if !found_reducer
                if negative
                    return false, changed
                end
                break
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
            if !is_monomial_reduction &&
                !GBElements.is_negative_disjoint(g, reducer, negative=negative)
                #@debug(
                #    "Found reducer but it is not negative disjoint",
                #    g, reducer, negative
                #)
                return true, changed
            end
            #Now apply the reduction and check if it is a zero reduction
            #@debug(
            #    "Reducing g by reducer", g, reducer, negative
            #)
            changed = true
            if has_order(gb) && !is_monomial_reduction
                reduced_to_zero = GBElements.reduce!(
                    g, reducer, order(gb), negative=negative
                )
            else
                reduced_to_zero = GBElements.reduce!(
                    g, reducer, negative=negative
                )
            end
            #@debug("Reduced g", result=g, reduced_to_zero)
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
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S}
) :: T where {T <: AbstractVector{Int}, S <: GBOrder}
    v = bs[i]
    w = bs[j]
    r = build(mem, v, w)
    orientate!(r, order(bs))
    return r
end

#It is critical to implement this efficiently, as it is called many times.
#isempty(intersect(b1, b2)) is not enough!
function disjoint(b1 :: BitSet, b2 :: BitSet)
    #Specialization of _matched_map of Julia / base / bitset.jl
    l1 = length(b1.bits)
    l2 = length(b2.bits)
    bdiff = b2.offset - b1.offset
    @inbounds for i = max(1, 1+bdiff):min(l1, l2+bdiff)
        if b1.bits[i] & b2.bits[i-bdiff] != 0
            return false
        end
    end
    return true
end

#disjoint(b1 :: BitSet, b2 :: BitSet) = isempty(intersect(b1, b2))

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

I also added new criteria for binary variables. In this case, a pair
(i, j) may be discarded if the positive binaries of i are disjoint from
the negative binaries of j, and vice-versa. This is because in these cases
the pair will generate some binomial with 2 or -2 in some coordinate.
These are never necessary in a Gröbner basis.
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    bs :: BinomialSet{T, S};
    use_binary_truncation :: Bool = true
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    #if is_minimization(bs)
    eliminate_by_criteria = !disjoint(negative_support(bs[i]), negative_support(bs[j])) ||
        disjoint(positive_support(bs[i]), positive_support(bs[j]))
    eliminate_by_truncation = use_binary_truncation && (
        !disjoint(positive_binaries(bs[i]), negative_binaries(bs[j])) ||
        !disjoint(negative_binaries(bs[i]), positive_binaries(bs[j])))
    return eliminate_by_criteria || eliminate_by_truncation
    #end
    #Maximization problem
    #return disjoint(negative_support(bs[i]), negative_support(bs[j])) ||
    #    !disjoint(positive_support(bs[i]), positive_support(bs[j]))
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
    is_truncated_groebner_basis(
    bs :: BinomialSet{T, S},
    truncate :: Function
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}

Returns true if and only if `bs` is a truncated Gröbner Basis with truncation
done using the `truncate` function.
"""
function is_truncated_groebner_basis(
    bs :: BinomialSet{T, S},
    truncate :: Function
) :: Bool where {T <: AbstractVector{Int}, S <: GBOrder}
    if isempty(bs)
        return true
    end
    mem = Vector{Int}(undef, length(bs[1]))
    for i in 1:length(bs)
        for j in 1:(i-1)
            s = BinomialPair(i, j)
            binomial = sbinomial(mem, s, bs)
            reduced_to_zero, _ = reduce!(binomial, bs)
            if !reduced_to_zero && !truncate(binomial)
                @info("Missing binomial in truncated GB", binomial,
                bs[i], bs[j], i, j)
                #Note we check truncation after reduction. This is often
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

function auto_reduce_once!(
    gb :: BinomialSet{T, S};
    current_index :: Int = 0
) :: Tuple{Int, Int} where {T <: AbstractVector{Int}, S <: GBOrder}
    removed = 0
    removed_before_index = 0
    i = length(gb)
    while i >= 1
        g = copy(gb[i]) #We reduce a copy of the element, otherwise, if it
        #is changed by reduction, it will stay in the reduction tree later
        reduced_to_zero, changed = reduce!(g, gb, skipbinomial=gb[i])
        if reduced_to_zero || changed
            if i < current_index
                removed_before_index += 1
                current_index -= 1
            end
            @debug("Autorreduced binomial", original=gb[i],
                result=(reduced_to_zero ? 0 : g),
                reduced_to_zero, changed
            )
            deleteat!(gb, i)
        end
        if reduced_to_zero
            removed += 1
        elseif changed
            #Elements have to be readded so that their S-pairs are
            #recomputed. Otherwise, the cancellation criterion may
            #break things.
            push!(gb, g)
        end
        i -= 1
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
        _, found_reducer = find_reducer(g, gb, reduction_tree(gb), skipbinomial=g)
        if found_reducer
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
    #TODO: Fix this later
    if has_signature(T)
        return
    end
    minimal_basis!(gb)
    if !is_minimization(gb)
        #Doesn't compute reduced basis for implicit representation
        #TODO: fix this in some later version!
        return
    end
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        while reducing
            h, found_reducer = find_reducer(
                g, gb, reduction_tree(gb), negative=true, skipbinomial=g
            )
            if found_reducer
                #The binomial has to be removed from the tree first, otherwise
                #its filter will already have been changed and it won't be
                #found in the tree.
                SupportTrees.remove_binomial!(gb.reduction_tree, g)
                GBElements.reduce!(g, h, order(gb), negative=true)
                SupportTrees.add_binomial!(gb.reduction_tree, g)
            else
                reducing = false
            end
        end
    end
end

end
