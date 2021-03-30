module Binomials

export Binomial, lattice_generator_binomial

using IPGBs.GBElements
using IPGBs.Orders

mutable struct Binomial <: GBElement
    element :: Vector{Int}
    cost :: Int
end

GBElements.cost(g :: Binomial) = g.cost
GBElements.fullform(g :: Binomial) = g.element

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

"""
Creates g - h in preallocated vector result.
"""
function GBElements.minus(
    result :: Vector{Int},
    g :: Binomial,
    h :: Binomial
) :: Binomial
    result .= g.element .- h.element
    new_cost = g.cost - h.cost
    return Binomial(result, new_cost)
end

function Base.copy(
    g :: Binomial
) :: Binomial
    return Binomial(copy(g.element), g.cost)
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
In addition to applying the reduction itself, we update the cost.

Returns true iff `binomial` reduced to 0.
"""
function GBElements.reduce!(
    g :: Binomial,
    h :: Binomial,
    order :: GBOrder;
    negative :: Bool = false
) :: Bool
    reduced_to_zero = GBElements.reduce!(g.element, h.element)
    if !reduced_to_zero #For efficiency, don't do the other stuff in this case
        g.cost -= h.cost
        orientate!(g, order)
    end
    return reduced_to_zero
end

function GBElements.opposite!(
    g :: Binomial
)
    g.element .= .-g.element
    g.cost = -g.cost
end

"""
Computes a Markov basis of `A` with `c` as cost matrix. This assumes the problem
is in the particular form given in Thomas and Weismantel (1997), Section 3.
"""
function lattice_generator_binomial(
    i :: Int,
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Array{Float64},
    u :: Vector{Int},
    order :: GBOrder;
    check_truncation :: Bool = true
) :: Union{Binomial, Nothing}
    #This assumes inequality + binary constraints
    #It is enough for my current experiments, but I should generalize this
    #The current matrix A has n + m rows, 2n + m cols. So n = cols - rows
    n = size(A, 2) - size(A, 1)
    m = size(A, 1) - n
    v = zeros(Int, n)
    v[i] = 1
    r = zeros(Int, n)
    r[i] = -1
    s = -copy(A[1:m, i])
    g = vcat(v, s, r)
    if ndims(c) == 1
        cost = c[i]
    else #c is two-dimensional
        @assert ndims(c) == 2
        cost = c[1, i]
    end
    generator = Binomial(g, cost)
    #The problem may be in minimization form, or have negative costs
    #Thus b may have negative cost, in that case we need to change its orientation
    orientate!(generator, order)
    #Check whether the Binomial should be truncated. If it should, just
    #return nothing instead
    if !check_truncation || isfeasible(generator, A, b, u)
        return generator
    end
    return nothing
end

end
