module Binomials

export initialize_binomials, Binomial, lattice_generator_binomial

using MIPMatrixTools.IPInstances

using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.Orders

element_end :: Int = 0
nonnegative_end :: Int = 0 #Index of the last non-negative variable of the current IP problem
bounded_end :: Int = 0 #Index of the last bounded variable of the current IP problem
cost_start :: Int = 0 #Index of the first cost value
cost_end :: Int = 0 #Index of the last cost value

function initialize_binomials(instance :: IPInstance, order :: MonomialOrder)
    global element_end = instance.n
    global nonnegative_end = instance.nonnegative_end
    global bounded_end = instance.bounded_end
    global cost_start = instance.n + 1
    global cost_end = instance.n + order.num_costs
end

struct Binomial <: GBElement
    #A binomial with costs. Costs are stored in the same vector for efficiency.
    #Data from indices cost_start to cost_end are the costs, data up to
    #element_end are the binomial's entries.
    data :: Vector{Int}
end

GBElements.data(g :: Binomial) = g.data
GBElements.element(g :: Binomial) = @views g.data[1:element_end]
GBElements.costs(g :: Binomial, :: GBOrder) = @views g.data[cost_start:cost_end]
GBElements.costs(g :: Binomial) = @views g.data[cost_start:cost_end]
GBElements.nonnegative(g :: Binomial) = @views g.data[1:nonnegative_end]
GBElements.bounded(g :: Binomial) = @views g.data[1:bounded_end]
GBElements.fullform(g :: Binomial) = @views g.data[1:element_end]

function Base.show(
    io :: IO,
    g :: Binomial
)
    print(io, element(g), " : c", costs(g))
end

function Base.:-(
    g :: Binomial,
    h :: Binomial
) :: Binomial
    new_element = g.data - h.data
    return Binomial(new_element)
end

"""
Creates g - h in preallocated vector result.
"""
function GBElements.minus(
    result :: Vector{Int},
    g :: Binomial,
    h :: Binomial
) :: Binomial
    result .= g.data .- h.data
    return Binomial(result)
end

function Base.copy(
    g :: Binomial
) :: Binomial
    return Binomial(copy(g.data))
end

#
# Implementation of the AbstractVector interface for Binomials
#

function Base.size(
    g :: Binomial
) :: Tuple
    return size(g.data)
end

function Base.getindex(
    g :: Binomial,
    i :: Int
) :: Int
    return g.data[i]
end

function Base.setindex!(
    g :: Binomial,
    v :: Int,
    i :: Int
)
    g.data[i] = v
end

function Base.length(
    g :: Binomial
) :: Int
    return length(g.data)
end

function Base.empty(g :: Binomial)
    return Binomial(empty(g.data))
end

#
# Other methods for Binomial
#

function GBElements.weight(g :: Binomial, w :: Vector{Float64})
    total_weight = 0.0
    for i in 1:nonnegative_end
        if g[i] > 0
            total_weight += g[i] * w[i]
        end
    end
    return total_weight
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
    if ndims(c) == 1
        cost = c[i]
    else #c is two-dimensional
        @assert ndims(c) == 2
        cost = c[1, i]
    end
    g = vcat(v, s, r, [cost])
    generator = Binomial(g)
    #The problem may be in minimization form, or have negative costs
    #Thus b may have negative cost, in that case we need to change its orientation
    orientate!(generator, order)
    #Check whether the Binomial should be truncated. If it should, just
    #return nothing instead
    if !check_truncation || simple_truncation(generator, A, b, u)
        return generator
    end
    return nothing
end

end
