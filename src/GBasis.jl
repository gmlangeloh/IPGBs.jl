module GBasis

using IPGBs.GBElements
using IPGBs.SupportTree

"""
This is specialized by MonomialOrder in case of Buchberger's algorithm and by
ModuleMonomialOrder in case of Signature-based algorithms.
"""
abstract type GBOrder end

struct MonomialOrder
    cost_matrix :: Array{Int, 2}
    #Probably should carry a tiebreaker around too
end

struct GrobnerBasis{T <: GBElement, S <: GBOrder} <: AbstractVector{T}
    basis :: Vector{T}
    order :: S
    reduction_tree :: SupportTree{T}
    matrix :: Array{Int, 2} #The data defining this toric ideal, matrix * x <= rhs
    rhs :: Vector{Int}
end

function Base.size(
    gb :: GrobnerBasis{T, S}
) :: Tuple where {T <: GBElement, S <: GBOrder}
    return size(gb.basis)
end

function Base.getindex(
    gb :: GrobnerBasis{T, S},
    i :: Int
) :: T where {T <: GBElement, S <: GBOrder}
    return gb.basis[i]
end

function Base.setindex(
    gb :: GrobnerBasis{T, S},
    g :: T,
    i :: Int
) where {T <: GBElement, S <: GBOrder}
    gb.basis[i] = g
end

function Base.length(
    gb :: GrobnerBasis{T, S}
) :: Int where {T <: GBElement, S <: GBOrder}
    return length(gb.basis)
end

function Base.push!(
    gb :: GrobnerBasis{T, S},
    g :: T
)
    push!(gb.basis, g)
    #TODO update this!
    push!(gb.module_ordering.generators, g)
    addbinomial!(gb.reduction_tree, gb[length(gb)])
end

end
