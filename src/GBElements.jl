"""
This module defines all the binomial data structures used in my implementations
of Buchberger's algorithm and Signature-based algorithms.
"""
module GBElements
export GBElement, Binomial, regular_reducible, degree_reducible, getfilter,
    lattice_generator_binomial, lt_tiebreaker, isfeasible, iszero

using IPGBs.FastBitSets

#
# Functions dealing with reducibility of any monomial or binomial-like vector-based
# structure
#

"""
Gets the filter of a binomial, that is, the list of indices of variables
appearing in its leading term.

If fullfilter = true, include the indices of variables appearing in its
trailing term.
"""
function getfilter(
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
    negative :: Bool = false
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
    #This is used in signature-based algorithms. Non-signature algorithms
    #simply return true here
    if !regular_reducible(reducer, g, gb)
        println("Singular reducer found")
        @show reducer reducer.signature
        @show g g.signature
        return false
    end
    return true
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

"""
Abstract type used for GB computations. It is meant to generalize both binomials
and binomials with signature, in order to simplify the implementation of
reduction algorithms.
"""
abstract type GBElement <: AbstractVector{Int} end

"""
Returns true iff `reducer` is regular wrt `g`, that is, sig(reducer) < sig(g).
This is used when checking for regular reductions in a signature-based
algorithm. For a non-signature-based algorithm, and thus a generic GBElement,
this is always true.
"""
function regular_reducible(
    reducer :: T,
    g :: T,
    gb :: S
) :: Bool where {T <: GBElement, S <: AbstractVector{T}}
    return true
end

#"""
#Checks whether element is smaller coordinate-wise than u only considering the
#given indices.
#"""
#function sparse_le(
#    indices :: Vector{Int},
#    element :: Vector{Int},
#    u :: Vector{Int}
#) :: Bool
#    return all(i -> abs(element[i]) <= u[i], indices)
#end

function le_upperbound(
    v :: T,
    u :: Vector{Int}
) :: Bool where {T <: GBElement}
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
    v :: T,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool where {T <: GBElement}
    head, tail = degrees(v, A)
    return le_coordinatewise(head, b) && le_coordinatewise(tail, b) &&
        le_upperbound(v, u)
        #sparse_le(v.head, v.element, u) && sparse_le(v.tail, v.element, u)
end

mutable struct Binomial <: GBElement
    element :: Vector{Int}
    cost :: Int
end

function Base.show(
    io :: IO,
    g :: Binomial
)
    print(io, g.element, " : c", g.cost)
end

function Base.:-(
    g :: Binomial,
    h :: Binomial
) :: Binomial
    new_element = g.element - h.element
    new_cost = g.cost - h.cost
    return Binomial(new_element, new_cost)
end

#
# Implementation of the AbstractVector interface for Binomials
#

function Base.size(
    g :: Binomial
) :: Tuple
    return size(g.element)
end

function Base.getindex(
    g :: Binomial,
    i :: Int
) :: Int
    return g.element[i]
end

function Base.setindex!(
    g :: Binomial,
    v :: Int,
    i :: Int
)
    g.element[i] = v
end

function Base.length(
    g :: Binomial
) :: Int
    return length(g.element)
end

#
# Other methods for Binomial
#

"""
This is only relevant when we consider the implicit representation of binomials.
For this explicit representation, we can always return true.
"""
function degree_reducible(
    g :: Binomial,
    h :: Binomial,
    negative :: Bool = false
) :: Bool
    return true
end

"""
In addition to applying the reduction itself, we update the cost.

Returns true iff `binomial` reduced to 0.
"""
function reduce!(
    g :: Binomial,
    h :: Binomial;
    negative :: Bool = false
) :: Bool
    reduced_to_zero = reduce!(g.element, h.element)
    #if reduced_to_zero
    #    return true
    #end
    #if g.cost > h.cost || (g.cost == h.cost && lt_tiebreaker(h, g))
    #    g.cost -= h.cost
    #else
    #    g.cost -= h.cost
    #    opposite!(g)
    #end
    g.cost -= h.cost
    orientate!(g)
    return reduced_to_zero
end

function opposite!(
    g :: Binomial
)
    g.element .= .-(g.element)
    g.cost = -g.cost
end

"""
Returns true iff g is oriented in a compatible way to grevlex with xn > xn-1 ...
> x1
"""
function grevlex(
    g :: Binomial
) :: Bool
    sum_g = sum(g)
    if sum_g > 0
        return true
    elseif sum_g < 0
        return false
    end
    i = 1
    while sum_g == 0 && i < length(g)
        sum_g -= g[i]
        if sum_g > 0
            return true
        elseif sum_g < 0
            return false
        end
        i += 1
    end
    return true
end

function orientate!(
    g :: Binomial
)
    #Applies tiebreaker by grevlex in case the cost is 0
    if g.cost < 0 || (g.cost == 0 && !grevlex(g))
        opposite!(g)
    end
end

function iszero(
    g :: Binomial
) :: Bool
    for gi in g
        if gi != 0
            return false
        end
    end
    return true
end

function degrees(
    g :: Binomial,
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
Computes a Markov basis of `A` with `c` as cost matrix. This assumes the problem
is in the particular form given in Thomas and Weismantel (1997), Section 3.
"""
function lattice_generator_binomial(
    i :: Int,
    A :: Array{Int, 2},
    c :: Array{Int}
) :: Binomial
    #This assumes inequality + binary constraints
    #It is enough for my current experiments, but I should generalize this
    m = size(A, 1)
    n = size(A, 2)
    v = zeros(Int, n)
    v[i] = 1
    r = zeros(Int, n)
    r[i] = -1
    s = -copy(A[:, i])
    g = vcat(v, s, r)
    if ndims(c) == 1
        cost = c[i]
    else #c is two-dimensional
        @assert ndims(c) == 2
        cost = c[1, i]
    end
    b = Binomial(g, cost)
    #The problem may be in minimization form, or have negative costs
    #Thus b may have negative cost, in that case we need to change its orientation
    orientate!(b)
    return b
end

"""
Computes bitsets with positive and negative supports of `g`.
"""
function supports(
    g :: Binomial
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

end
