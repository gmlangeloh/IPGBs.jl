"""
This module defines all the binomial data structures used in my implementations
of Buchberger's algorithm and Signature-based algorithms.
"""
module GBElements
export GBElement, regular_reducible, iszero, degree_reducible, getfilter

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

TODO I need something else besides fullfilter to tell if I should use degree_red
and reg_red, because my new Binomial structure shouldn't use fullfilter
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
        #This is used in signature-based algorithms. Non-signature algorithms
        #simply return true here
        if !regular_reducible(reducer, g, gb)
            println("Singular reducer found")
            @show reducer reducer.signature
            @show g g.signature
            return false
        end
    else
        for i in filter
            if sign * g[i] < reducer[i]
                return false
            end
        end
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
    while reducer[i] <= 0
        i += 1
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
    g.cost -= h.cost
    return reduced_to_zero
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

function lattice_generator(
    i :: Int,
    A :: Array{Int, 2},
    c :: Array{Int}
) :: Binomial
    #TODO implement this...
end

end
