"""
This module defines all the binomial data structures used in my implementations
of Buchberger's algorithm and Signature-based algorithms.

TODO make GBElements a consistent interface
"""
module GBElements
export GBElement, degree_reducible, filter, lt_tiebreaker, isfeasible, is_zero, leading_term, head, has_signature, singular_top_reducible, signature_reducible

using IPGBs.FastBitSets

"""
Abstract type used for GB computations. It is meant to generalize both binomials
and binomials with signature, in order to simplify the implementation of
reduction algorithms.
"""
abstract type GBElement <: AbstractVector{Int} end

#
# Hard contract: a GBElement must implement at least the following functions.
#

#TODO will I even have anything here?

#
# Soft contract: concrete GBElements may reimplement these if necessary
#

has_signature(g :: AbstractVector{Int}) = false
has_signature(g :: GBElement) = false

"""
Computes the leading term of this GBElement as a vector.
"""
function leading_term(
    g :: GBElement
) :: Vector{Int}
    lt = zeros(Int, length(g))
    for i in 1:length(g)
        if g[i] > 0
            lt[i] = g[i]
        end
    end
    return lt
end

"""
The indices of the positive support of `g`. Or indices of the support of
leading_term(g)

TODO could probably turn this into an iterator instead, would be more efficient
"""
function head(
    g :: GBElement
) :: Vector{Int}
    head = Int[]
    for i in 1:length(g)
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
    g :: GBElement
) :: Bool
    for gi in g
        if gi != 0
            return false
        end
    end
    return true
end

function opposite!(
    g :: GBElement
)
    for i in 1:length(g)
        g[i] = -g[i]
    end
end

function degrees(
    g :: GBElement,
    A :: Array{Int, 2}
) :: Tuple{Vector{Int}, Vector{Int}}
    #Get the positive and negative parts of g
    positive_g = Int[]
    negative_g = Int[]
    for i in 1:length(g)
        if g[i] > 0
            push!(positive_g, g[i])
            push!(negative_g, 0)
        else
            push!(negative_g, g[i])
            push!(positive_g, 0)
        end
    end
    #Compute the degrees of positive_g and negative_g by A
    return A * positive_g, A * negative_g
end

"""
Computes bitsets with positive and negative supports of `g`.
"""
function supports(
    g :: GBElement
) :: Tuple{FastBitSet, FastBitSet}
    pos_supp = Int[]
    neg_supp = Int[]
    for i in 1:length(g)
        if g[i] > 0
            push!(pos_supp, i)
        elseif g[i] < 0
            push!(neg_supp, i)
        end
    end
    bitset_length = length(g)
    return makebitset(bitset_length, pos_supp), makebitset(bitset_length, neg_supp)
end

"""
This is only relevant when we consider the implicit representation of binomials.
For this explicit representation, we can always return true.
"""
function degree_reducible(
    g :: T,
    h :: T;
    negative :: Bool = false
) :: Bool where {T <: GBElement}
    return true
end

#For generic GBElements which don't has_signature(g), these are preset.
signature_reducible(g :: GBElement, reducer_sig, gb) = true
singular_top_reducible(g :: GBElement, reducer_sig) = false

#
# Additional GBElement logic
#

"""
Checks whether v is bounded coordinate by coordinate by u.
"""
function le_upperbound(
    v :: GBElement,
    u :: Vector{Int}
) :: Bool
    for i in 1:length(v)
        if v[i] > 0 && v[i] > u[i]
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
    return all(i -> v1[i] <= v2[i], keys(v1))
end

"""
Returns true iff v should be considered for reduction in a truncated GB
algorithm.
"""
function isfeasible(
    v :: GBElement,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool
    head, tail = degrees(v, A)
    return le_coordinatewise(head, b) && le_coordinatewise(tail, b) &&
        le_upperbound(v, u)
        #sparse_le(v.head, v.element, u) && sparse_le(v.tail, v.element, u)
end

#
# Functions dealing with reducibility of any monomial or binomial-like vector-based
# structure. Work for any GBElement or, more generically, any AbstractVector{Int}
# This means they support monomials as well
#

"""
Gets the filter of a binomial, that is, the list of indices of variables
appearing in its leading term.

If fullfilter = true, include the indices of variables appearing in its
trailing term.
"""
function filter(
    binomial :: T;
    fullfilter :: Bool = false
) :: Vector{Int} where {T <: AbstractVector{Int}}
    filter = Int[]
    for i in 1:length(binomial)
        if binomial[i] > 0 || (fullfilter && binomial[i] != 0)
            push!(filter, i)
        end
    end
    return filter
end

"""
Checks whether `reducer` divides `g`, using the filter of `reducer` for
efficiency. When fullfilter = true, checks if g.head >= reducer.head and
g.tail >= reducer.tail coordinate-wise, while also checking a degree criterion.
"""
function reduces(
    g :: T,
    filter :: Vector{Int},
    reducer :: T,
    gb :: S;
    fullfilter :: Bool = true,
    negative :: Bool = false,
    params :: Dict = Dict()
) :: Bool where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
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
    if has_signature(g)
        #1. Compute the reduction factor, for efficiency. It will be used in all
        #of the following tests
        quotient = monomial_quotient(g, reducer)
        reducer_sig = quotient * reducer.signature
        #2. Check if singular. If so, we need to stop looking for reducers
        if singular_top_reducible(g, reducer_sig)
            #We set this to true to signal the support tree to stop looking for
            #reducers
            params["is_singular"] = true
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

function monomial_quotient(
    binomial :: T,
    reducer :: T
) :: Vector{Int} where {T <: GBElement}
    n = length(reducer)
    quotient = zeros(Int, n)
    for i in head(binomial)
        quotient[i] = binomial[i] - reducer[i]
        #We do assume that reducer divides binomial
        @assert quotient[i] >= 0
    end
    return quotient
end

"""
Finds the maximum k such that k * reducer <= binomial coordinate-wise.
"""
function reduction_factor(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Int where {T <: AbstractVector{Int}}
    i = 1
    while i <= length(reducer) && reducer[i] <= 0
        i += 1
    end
    if i > length(reducer)
        return negative ? -1 : 1
    end
    factor = Int(floor(binomial[i] / reducer[i]))
    if (!negative && factor == 1) || (negative && factor == -1)
        return factor
    end
    i += 1
    while i <= length(reducer)
        if reducer[i] > 0
            newfactor = Int(floor(binomial[i] / reducer[i]))
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
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}}
    factor = reduction_factor(binomial, reducer, negative=negative)
    reduced_to_zero = true
    if !negative && factor == 1
        for i in 1:length(binomial)
            binomial[i] -= reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    elseif negative && factor == -1
        for i in 1:length(binomial)
            binomial[i] += reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    else
        for i in 1:length(binomial)
            binomial[i] -= factor * reducer[i]
            if binomial[i] != 0
                reduced_to_zero = false
            end
        end
    end
    return reduced_to_zero
end

"""
Returns true iff g is smaller than h in the tiebreaking order (grevlex).
Assumes that g.cost == h.cost, that is, that there is a tie.

I don't really understand why this works, it doesn't look like grevlex to me.
But it gives precisely the same results as 4ti2, so I guess I'll keep it.
Commented below is an implementation which does not give the same results as
4ti2, but makes more sense to me.
"""
function lt_tiebreaker(
    g :: T,
    h :: T
) :: Bool where {T <: AbstractVector{Int}}
    @assert g.cost == h.cost
    gsmaller :: Int = 0 #-1 when g < h, 0 when g = h, 1 when g > h
    sum_g :: Int = 0
    sum_h :: Int = 0
    # Compute cumulative sums for g.element and h.element
    #the smallest one is the one with lowest cumulative sum at the farthest
    #point where they don't tie.
    for i in 1:length(g)
        sum_g += g[i]
        sum_h += h[i]
        if sum_g < sum_h
            gsmaller = -1
            break
        elseif sum_g > sum_h
            gsmaller = 1
            break
        end
    end
    return gsmaller == -1 ? true : false
    #for i in 1:length(g)
    #    sum_g += g[i]
    #    sum_h += h[i]
    #end
    #if sum_g < sum_h
    #    return true
    #elseif sum_g > sum_h
    #    return false
    #end
    #for i in 1:length(g)
    #    sum_g -= g[i]
    #    sum_h -= h[i]
    #    if sum_g < sum_h
    #        return true
    #    elseif sum_g > sum_h
    #        return false
    #    end
    #end
    #return false #If they are equal wrt grevlex at this point, g == h
end

end
