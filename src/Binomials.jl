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
binary_variables :: Vector{Bool} = Bool[] #1 if the variable is binary in the IP

const EMPTY_FILTER :: Vector{Int} = Int[]
const EMPTY_BITSET :: FastBitSet = FastBitSet(0)

function initialize_binomials(instance :: IPInstance, order :: MonomialOrder)
    global element_end = instance.n
    global nonnegative_end = instance.nonnegative_end
    global bounded_end = instance.bounded_end
    global cost_start = instance.n + 1
    global cost_end = instance.n + order.num_costs
    global binary_variables = instance.binaries
end

mutable struct Binomial <: GBElement
    #A binomial with costs. Costs are stored in the same vector for efficiency.
    #Data from indices cost_start to cost_end are the costs, data up to
    #element_end are the binomial's entries.
    data :: Vector{Int}

    #The following fields are used to efficiently obtain the indices
    #of non-zero entries of the binomial.
    computed_supports :: Bool
    positive_support :: FastBitSet
    negative_support :: FastBitSet
    positive_filter :: Vector{Int}
    negative_filter :: Vector{Int}

    positive_binaries :: FastBitSet
    negative_binaries :: FastBitSet

    function Binomial(v :: Vector{Int})
        #By default, we do not compute the supports. They are computed lazily.
        #Only binomials added to a GB need them.
        pos_filter = EMPTY_FILTER
        neg_filter = EMPTY_FILTER
        pos_supp = EMPTY_BITSET
        neg_supp = EMPTY_BITSET
        pos_binaries = EMPTY_BITSET
        neg_binaries = EMPTY_BITSET
        new(v, false, pos_supp, neg_supp, pos_filter, neg_filter,
            pos_binaries, neg_binaries
        )
    end
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

function GBElements.compute_supports(g :: Binomial)
    if g.computed_supports
        return
    end
    pos_supp, neg_supp, pos_filter, neg_filter = GBElements.supports(g)
    g.positive_support = pos_supp
    g.negative_support = neg_supp
    g.positive_filter = pos_filter
    g.negative_filter = neg_filter
    g.computed_supports = true
end

function GBElements.compute_binaries(g :: Binomial)
    pos_bin = Int[]
    neg_bin = Int[]
    for i in eachindex(element(g))
        if binary_variables[i]
            if g[i] > 0
                push!(pos_bin, i)
            elseif g[i] < 0
                push!(neg_bin, i)
            end
        end
    end
    g.positive_binaries = FastBitSet(length(element(g)), pos_bin)
    g.negative_binaries = FastBitSet(length(element(g)), neg_bin)
end

GBElements.positive_support(g :: Binomial) = g.positive_support
GBElements.negative_support(g :: Binomial) = g.negative_support
GBElements.positive_filter(g :: Binomial) = g.positive_filter
GBElements.negative_filter(g :: Binomial) = g.negative_filter
GBElements.positive_binaries(g :: Binomial) = g.positive_binaries
GBElements.negative_binaries(g :: Binomial) = g.negative_binaries

end
