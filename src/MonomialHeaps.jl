module MonomialHeaps

using DataStructures
using IPGBs.Orders

export WeightedMonomial, MonomialHeap

"""
Monomial with `weight` given by a monomial order, used for fast comparison in
a priority queue or heap. It also includes a `count` to keep track of how many
times this monomial has been generated. FGLM can use this to check whether
some monomial is underneath a GB staircase efficiently.
"""
mutable struct WeightedMonomial
    monomial :: Vector{Int}
    weight :: Float64
    count :: Int
    divisor_nf :: Vector{Int} #The normal form of a monomial dividing `monomial` when available
    var_index :: Int #A variable such that divisor[var_index] + 1 == monomial[var_index] and divisor[i] == monomial[i] for all other i. Set to 0 in case this isn't possible
    support_size :: Int
end

function WeightedMonomial(
    m::Vector{Int},
    o::S,
    divisor::Vector{Int},
    i::Int,
    supp_size :: Int
) where {S<:GBOrder}
    return WeightedMonomial(m, order_cost(o, m), 1, divisor, i, supp_size)
end

Base.isless(m1::WeightedMonomial, m2::WeightedMonomial) = m1.weight < m2.weight

"""
Keeps a priority list / heap of monomials, ordered by weight, without
repetitions. The number of times any monomial is inserted is stored in the
WeightedMonomial structure.
"""
struct MonomialHeap{T <: GBOrder}
    heap :: BinaryMinHeap{WeightedMonomial}
    set :: Dict{Vector{Int}, WeightedMonomial}
    order :: T

    function MonomialHeap{T}(o :: T) where {T <: GBOrder}
        h = BinaryMinHeap{WeightedMonomial}()
        s = Dict{Vector{Int}, WeightedMonomial}()
        new(h, s, o)
    end
end

function MonomialHeap(
    o :: T,
    monoms :: Vector{Vector{Int}}
) where {T <: GBOrder}
    heap = MonomialHeap{T}(o)
    for m in monoms
        push!(heap, m)
    end
    return heap
end

Base.isempty(h :: MonomialHeap) = isempty(h.heap)

function Base.push!(
    h :: MonomialHeap,
    m :: Vector{Int},
    divisor :: Vector{Int},
    var_index :: Int,
    supp_size :: Int
)
    if haskey(h.set, m)
        wm = h.set[m]
        wm.count += 1
    else
        wm = WeightedMonomial(m, h.order, divisor, var_index, supp_size)
        push!(h.heap, wm)
        h.set[m] = wm
    end
end

Base.push!(h::MonomialHeap, m::Vector{Int}) =
    push!(h,m,zeros(Int,length(m)),0, 0)

function Base.pop!(h :: MonomialHeap)
    wm = pop!(h.heap)
    delete!(h.set, wm.monomial)
    return wm
end

end
