""" This module defines all the binomial data structures used in my implementations
of Buchberger's algorithm and Signature-based algorithms.

TODO: make GBElements a consistent interface
"""
module GBElements
#TODO: this is way too long, clean it up or at least break it into more exports
export GBElement, degree_reducible, filter, simple_truncation, is_zero, leading_term, head, has_signature, singular_top_reducible, signature_reducible, fullform, cost, CriticalPair, BinomialPair, first, second, build, is_implicit, orientate!, is_negative_disjoint, model_truncation, truncate, ipgbs_form, to_gbelement, weight, data, element, costs, bounded, nonnegative

using IPGBs.FastBitSets
using IPGBs.Orders
using IPGBs.SolverTools

using JuMP

"""
Abstract type used for GB computations. It is meant to generalize both binomials
and binomials with signature, in order to simplify the implementation of
reduction algorithms.
"""
abstract type GBElement <: AbstractVector{Int} end

#
# Hard contract: a GBElement must implement at least the following functions.
#

fullform(:: GBElement) :: Vector{Int} = error("Not implemented.")
minus(:: Vector{Int}, :: GBElement, :: GBElement) = error("Not implemented.")
weight(:: GBElement, :: Vector{Float64}) = error("Not implemented")

#
# Soft contract: concrete GBElements may reimplement these if necessary
#

has_signature(:: Type{<: AbstractVector{Int}}) = false
is_implicit(:: Type{<: AbstractVector{Int}}) = false

data(v :: AbstractVector{Int}) = v
element(v :: AbstractVector{Int}) = v
nonnegative(v :: AbstractVector{Int}) = v
bounded(v :: AbstractVector{Int}) = v
costs(v :: AbstractVector{Int}, o :: GBOrder) = order_costs(o, v)

"""
Turns a vector `v` into a GBElement of type `S`. Currently, Binomials are the
only subtype supported.
"""
function to_gbelement(
    v :: Vector{Int},
    order :: T,
    S :: DataType,
) where {T <: GBOrder}
    costs = round.(Int, order_costs(order, v))
    b = S([v; costs])
    orientate!(b, order)
    return b
end

"""
Turns a 4ti2 GB into my GB format as a vector of integer vectors.
"""
function ipgbs_form(
    fourti2_gb :: Matrix{Int}
) :: Vector{Vector{Int}}
    m = size(fourti2_gb, 1)
    gb = Vector{Int}[]
    for i in 1:m
        gb_elem = fourti2_gb[i,:]
        push!(gb, gb_elem)
    end
    return gb
end

"""
Computes the leading term of this GBElement as a vector.
"""
function leading_term(
    g :: T
) :: Vector{Int} where {T <: AbstractVector{Int}}
    lt = zeros(Int, length(element(g)))
    for i in eachindex(element(g))
        if g[i] > 0
            lt[i] = g[i]
        end
    end
    return lt
end

"""
The indices of the positive support of `g`. Or indices of the support of
leading_term(g)

TODO: could probably turn this into an iterator instead, would be more efficient
"""
function head(
    g :: GBElement
) :: Vector{Int}
    head = Int[]
    for i in eachindex(element(g))
        if g[i] > 0
            push!(head, i)
        end
    end
    return head
end

"""
Returns true iff this binomial is zero, that is, all of its coordinates are 0.
"""
function is_zero(
    g :: T
) :: Bool where {T <: AbstractVector{Int}}
    for gi in g
        if gi != 0
            return false
        end
    end
    return true
end

opposite!(g :: T) where {T <: AbstractVector{Int}} = g .= .-g

function orientate!(
    g :: T,
    order :: GBOrder
) where {T <: AbstractVector{Int}}
    cs = costs(g, order)
    if is_inverted(order, element(g), cs)
        GBElements.opposite!(g)
    end
end

function degrees(
    g :: T,
    A :: Array{Int, 2}
) :: Tuple{Vector{Int}, Vector{Int}} where {T <: AbstractVector{Int}}
    #Get the positive and negative parts of g
    n = length(element(g))
    positive_g = zeros(Int, n)
    negative_g = zeros(Int, n)
    for i in 1:n
        if g[i] > 0
            positive_g[i] = g[i]
        else
            negative_g[i] = g[i]
        end
    end
    #Compute the degrees of positive_g and negative_g by A
    return A * positive_g, A * negative_g
end

"""
Computes bitsets with positive and negative supports of `g`.
"""
function supports(
    g :: AbstractVector{Int}
) :: Tuple{FastBitSet, FastBitSet}
    pos_supp = Int[]
    neg_supp = Int[]
    for i in eachindex(nonnegative(g))
        if g[i] > 0
            push!(pos_supp, i)
        end
    end
    for i in eachindex(bounded(g))
        if g[i] < 0
            push!(neg_supp, i)
        end
    end
    pos_bitset_length = length(nonnegative(g))
    neg_bitset_length = length(bounded(g))
    return FastBitSet(pos_bitset_length, pos_supp), FastBitSet(neg_bitset_length, neg_supp)
end

"""
Computes bitsets with positive and negative supports of every element in `B`.
"""
function supports(
    B :: Vector{T}
) :: Tuple{Vector{FastBitSet}, Vector{FastBitSet}} where {T <: AbstractVector{Int}}
    pos_supps = FastBitSet[]
    neg_supps = FastBitSet[]
    for g in B
        p, n = supports(g)
        push!(pos_supps, p)
        push!(neg_supps, n)
    end
    return pos_supps, neg_supps
end

"""
This is only relevant when we consider the implicit representation of binomials.
For this explicit representation, we can always return true.
"""
function degree_reducible(
    g :: T,
    h :: T;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    return true
end

#For generic GBElements which don't has_signature(g), these are preset.
signature_reducible(g :: GBElement, reducer_sig, gb) = true
singular_top_reducible(g :: GBElement, reducer_sig) = false

#
# Additional GBElement logic
#

"""
Returns true iff `binomial` should be truncated according to the given
`truncation_type`.
"""
function truncate(
    binomial :: T,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Union{Int, Nothing}},
    model :: JuMP.Model,
    model_constraints :: Vector{JuMP.ConstraintRef},
    should_truncate :: Bool,
    truncation_type :: Symbol
) :: Bool where {T <: AbstractVector{Int}}
    if !should_truncate
        return false
    end
    if truncation_type == :Simple
        truncated = !simple_truncation(binomial, A, b, u)
    else #algorithm.truncation_type == :Model
        truncated = !model_truncation(
            binomial, A, b, model, model_constraints
        )
    end
    return truncated
end

"""
Checks whether the trailing terms of `g` and `h` are disjoint. In case
negative = true, checks whether the leading term of `g` is disjoint
with the trailing term of `h`.
"""
function is_negative_disjoint(
    g :: T,
    h :: T;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    sign = negative ? -1 : 1
    for i in eachindex(bounded(g))
        if sign * g[i] < 0 && h[i] < 0
            return false
        end
    end
    return true
end

"""
Checks whether v is bounded coordinate by coordinate by u.
"""
function le_upperbound(
    v :: T,
    u :: Vector{Union{Int, Nothing}}
) :: Bool where {T <: AbstractVector{Int}}
    for i in eachindex(element(v))
        if isnothing(u[i])
            continue
        elseif v[i] > 0 && v[i] > u[i]
            return false
        elseif v[i] < 0 && -v[i] > u[i]
            return false
        end
    end
    return true
end

"""
Returns true iff v1 <= v2 coordinate-wise.
"""
function le_coordinatewise(
    v1 :: Vector{Int},
    v2 :: Vector{Int}
) :: Bool
    return all(i -> v1[i] <= v2[i], keys(element(v1)))
end

"""
Returns true iff v should be considered for reduction in a truncated GB
algorithm.
"""
function simple_truncation(
    v :: T,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Union{Int, Nothing}}
) :: Bool where {T <: AbstractVector{Int}}
    head, tail = degrees(v, A)
    return le_coordinatewise(head, b) && le_coordinatewise(tail, b) &&
        le_upperbound(v, u)
end

"""
Checks feasibility of:

Ax = b - Av
0 <= x <= u

with x either a vector of real variables or integer variables, whatever is
defined in `model`. In order to do this efficiently, updates the RHS of the
previously built `model` by modifying `constraints`.

Returns true iff `v` is feasible for the above model, which means it should NOT
be truncated.

TODO: LP truncation is slightly slower than I hoped. IP truncation is very slow.
Can I implement the model more efficiently somehow?
"""
function model_truncation(
    v :: T,
    A :: Array{Int, 2},
    b :: Vector{Int},
    model :: JuMP.Model,
    constraints :: Vector{ConstraintRef}
) :: Bool where  {T <: AbstractVector{Int}}
    SolverTools.update_feasibility_model_rhs(constraints, A, b, leading_term(v))
    return SolverTools.is_feasible(model)
end

#
# Functions dealing with reducibility of any monomial or binomial-like vector-based
# structure. Work for any GBElement or, more generically, any AbstractVector{Int}
# This means they support monomials as well
#

"""
    filter(
    binomial :: T;
    fullfilter :: Bool = false
) :: Vector{Int} where {T <: AbstractVector{Int}}

The filter of a binomial, that is, the list of indices of variables
appearing in its leading term.

If fullfilter = true, include the indices of variables appearing in its
trailing term.
"""
function filter(
    binomial :: T;
    fullfilter :: Bool = false
) :: Vector{Int} where {T <: AbstractVector{Int}}
    filter = Int[]
    for i in eachindex(nonnegative(binomial))
        if binomial[i] > 0 || (fullfilter && binomial[i] != 0)
            push!(filter, i)
        end
    end
    return filter
end

"""
    reduces(
    g :: P,
    filter :: Vector{Int},
    reducer :: T,
    gb :: S;
    fullfilter :: Bool = true,
    negative :: Bool = false,
    is_singular :: Ref{Bool} = Ref(false)
) :: Bool where {P <: AbstractVector{Int}, T <: AbstractVector{Int}, S <: AbstractVector{T}}

Checks whether `reducer` divides `g`, using the filter of `reducer` for
efficiency. When fullfilter = true, checks if g.head >= reducer.head and
g.tail >= reducer.tail coordinate-wise, while also checking a degree criterion.
"""
function reduces(
    g :: P,
    filter :: Vector{Int},
    reducer :: T,
    gb :: S;
    fullfilter :: Bool = true,
    negative :: Bool = false
) :: Bool where {P <: AbstractVector{Int}, T <: AbstractVector{Int}, S <: AbstractVector{T}}
    sign :: Int = negative ? -1 : 1
    if fullfilter
        for i in filter
            if sign * g[i] * reducer[i] < 0 #Different signs, doesn't divide
                return false
            elseif abs(g[i]) < abs(reducer[i]) #reducer doesn't divide
                return false
            end
        end
        #Checks the truncation criterion of Thomas and Weismantel
        #This is only relevant because certain variables are implicit in this
        #representation and so this check is necessary to guarantee these variables
        #are compatible for division
        if !degree_reducible(reducer, g, negative=negative)
            return false
        end
    else
        for i in filter
            if sign * g[i] < reducer[i]
                return false
            end
        end
    end
    #This is used in signature-based algorithms. Non-signature algorithms just
    #skip this part
    if has_signature(P)
        #1. Compute the reduction factor, for efficiency. It will be used in all
        #of the following tests
        quotient = monomial_quotient(g, reducer)
        reducer_sig = quotient * reducer.signature
        #2. Check if singular. If so, we need to stop looking for reducers
        if singular_top_reducible(g, reducer_sig)
            #We set this to true to signal the support tree to stop looking for
            #reducers
            #is_singular[] = true
            return true #This stops the search. The current reducer won't matter
        end
        #3. Check if signature-reducible. This means reducer_sig < sig(g)
        if !signature_reducible(g, reducer_sig, gb)
            #In this case, we cannot s-reduce g by reducer, but we should keep
            #looking for reducers of g
            return false
        end
    end
    return true
end

"""
    monomial_quotient(
    binomial :: T,
    reducer :: T
) :: Vector{Int} where {T <: GBElement}

TBW
"""
function monomial_quotient(
    binomial :: T,
    reducer :: T
) :: Vector{Int} where {T <: GBElement}
    n = length(element(reducer))
    quotient = zeros(Int, n)
    for i in head(binomial)
        quotient[i] = binomial[i] - reducer[i]
        #We do assume that reducer divides binomial
        @assert quotient[i] >= 0
    end
    return quotient
end

"""
    reduction_factor(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Int where {T <: AbstractVector{Int}}

Finds the maximum k such that k * reducer <= binomial coordinate-wise.
"""
function reduction_factor(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Int where {T <: AbstractVector{Int}}
    i = 1
    n = length(element(reducer))
    while i <= n && reducer[i] <= 0
        i += 1
    end
    if i > n
        return negative ? -1 : 1
    end
    factor = Int(trunc(binomial[i] / reducer[i]))
    if (!negative && factor == 1) || (negative && factor == -1)
        return factor
    end
    i += 1
    while i <= n
        if reducer[i] > 0
            newfactor = Int(trunc(binomial[i] / reducer[i]))
            if (!negative && newfactor < factor) ||
                (negative && newfactor > factor)
                factor = newfactor
                if (!negative && factor == 1) || (negative && factor == -1)
                    return factor
                end
            end
        end
        i += 1
    end
    return factor
end

"""
Reduce `binomial` by `reducer` as many times as possible. Assumes `reducer`
reduces `binomial` at least once.

Returns true iff `binomial` reduced to 0.
"""
function reduce!(
    b :: T,
    r :: T,
    o :: GBOrder;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    reduced_to_zero = reduce!(b, r, negative=negative)
    if !reduced_to_zero
        orientate!(b, o)
    end
    return reduced_to_zero
end

#TODO: Optimize! Why is there no vectorization here? This is a mess!
function reduce!(
    binomial :: T,
    reducer :: T,
    factor :: Int;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    reduced_to_zero = true
    if !negative && factor == 1
        for i in eachindex(binomial)
            binomial[i] -= reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    elseif negative && factor == -1
        for i in eachindex(binomial)
            binomial[i] += reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    else
        for i in eachindex(binomial)
            binomial[i] -= factor * reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    end
    return reduced_to_zero
end

function reduce!(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    factor = reduction_factor(binomial, reducer, negative=negative)
    return reduce!(binomial, reducer, factor, negative=negative)
end

#
# Implementation of Critical Pairs. Can be extended by Signature algorithms.
#

abstract type CriticalPair end

#CriticalPair interface
first(:: CriticalPair) :: Int = error("Not implemented.")
second(:: CriticalPair) :: Int = error("Not implemented.")

struct BinomialPair <: CriticalPair
    i :: Int
    j :: Int
end

first(pair :: BinomialPair) = pair.i
second(pair :: BinomialPair) = pair.j

function build(
    mem :: Vector{Int},
    u :: T,
    v :: T,
    pair :: CriticalPair #This is used for more complicated CriticalPairs
) :: T where {T <: AbstractVector{Int}}
    return minus(mem, u, v)
end

end
