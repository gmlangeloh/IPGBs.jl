module Binomials

export Binomial, lattice_generator_binomial

using IPGBs.GBElements

mutable struct Binomial <: GBElement
    element :: Vector{Int}
    cost :: Int
end

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
    h :: Binomial;
    negative :: Bool = false
) :: Bool
    reduced_to_zero = GBElements.reduce!(g.element, h.element)
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

function GBElements.opposite!(
    g :: Binomial
)
    g.element .= .-g.element
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
        GBElements.opposite!(g)
    end
end

"""
Computes a Markov basis of `A` with `c` as cost matrix. This assumes the problem
is in the particular form given in Thomas and Weismantel (1997), Section 3.
"""
function lattice_generator_binomial(
    i :: Int,
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Array{Int},
    u :: Vector{Int};
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
    orientate!(generator)
    #Check whether the Binomial should be truncated. If it should, just
    #return nothing instead
    if !check_truncation || isfeasible(generator, A, b, u)
        return generator
    end
    return nothing
end

end
