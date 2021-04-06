module Orders

export GBOrder, MonomialOrder, is_inverted

using IPGBs.GBTools

"""
This is specialized by MonomialOrder in case of Buchberger's algorithm and by
ModuleMonomialOrder in case of Signature-based algorithms.
"""
abstract type GBOrder <: Base.Ordering end

is_inverted(:: GBOrder, :: AbstractVector{Int}) = error("Not implemented.")

"""
Returns a new matrix with the rows of C plus new rows consisting of
a grevlex tiebreaker, if C might not specify a full monomial order.
"""
function build_order(
    C :: Array{Float64, 2},
    A :: Array{Int, 2},
    b :: Vector{Int}
) :: Array{Float64, 2}
    m = size(C, 1)
    n = size(C, 2)
    if m < n #Insufficient tiebreaking in C
        tiebreaker = GBTools.revlex_matrix(n)
        full_matrix = vcat(C, tiebreaker)
    else
        full_matrix = C
    end
    #Make the first row strictly positive
    if any(full_matrix[1, j] < 0 for j in 1:n)
        d = GBTools.positive_row_span(A, b)
        lambda = 1 + maximum(-full_matrix[1, j] / d[j] for j in 1:n)
        for j in 1:n
            full_matrix[1, j] += lambda * d[j]
        end
    end
    #Store the transpose to exploit the fact Julia stores matrices
    #in column-major form. This helps performance of cmp and lt.
    return full_matrix'
end

"""
Implements a monomial order from a given cost matrix, including a grevlex
tiebreaker if necessary.
"""
mutable struct MonomialOrder <: GBOrder
    cost_matrix :: Array{Float64, 2}
    tiebreaker :: Symbol
    #Probably should carry a tiebreaker around too
    MonomialOrder(costs :: Array{Float64, 2}, A, b) = new(build_order(costs, A, b), :invlex)
end

Base.length(o :: MonomialOrder) = size(o.cost_matrix, 2)

#This is implemented this way for locality and performance. It is only used
#in cmp / lt below regardless.
Base.getindex(o :: MonomialOrder, i, j) = o.cost_matrix[j, i]

"""
Efficient comparison of vectors with respect to a monomial order.
"""
function Base.cmp(
    order :: MonomialOrder,
    v1 :: T,
    v2 :: T
) :: Int where {T <: AbstractVector{Int}}
    @assert length(v1) == length(v2)
    for i in 1:length(order)
        #Compute o' * (v1 - v2) without allocations
        s = 0
        for j in 1:length(v1)
            s += order[i, j] * (v1[j] - v2[j])
        end
        if s < 0
            return -1
        elseif s > 0
            return 1
        end #s == 0, tie. Try next vector in order
    end
    return 0 #tied completely, v1 == v2.
end

function Base.lt(
    order :: MonomialOrder,
    v1 :: T,
    v2 :: T
) where {T <: AbstractVector{Int}}
    return cmp(order, v1, v2) == -1
end

"""
Returns true iff the trailing and leading terms of `v` are inverted with
respect to `order`.
"""
function is_inverted_generic(
    order :: MonomialOrder,
    v :: T
) :: Bool where {T <: AbstractVector{Int}}
    for i in 1:length(order)
        s = 0
        for j in 1:length(v)
            s += order[i, j] * v[j]
        end
        if s < 0
            return true
        elseif s > 0
            return false
        end
    end
    return false #This really should not happen
end

"""
`v` is inverted according to invlex iff its first non-zero entry is positive.
"""
function is_inverted_invlex(
    v :: T
) :: Bool where {T <: AbstractVector{Int}}
    i = 1
    while i <= length(v) && v[i] == 0
        i += 1
    end
    if i > length(v)
        return false
    end
    return v[i] > 0
end

"""
Efficiently returns whether v's leading and trailing terms are inverted.
"""
function is_inverted(
    order :: MonomialOrder,
    v :: T,
    cost :: Int
) :: Bool where {T <: AbstractVector{Int}}
    if cost < 0
        return true
    elseif cost > 0
        return false
    elseif order.tiebreaker == :invlex
        return is_inverted_invlex(v)
    end
    return is_inverted_generic(order, v)
end

"""
Changes the order stored in MonomialOrder to a new order. Rebuilds
tiebreaker if necessary.
"""
function change_ordering!(
    order :: MonomialOrder,
    new_order :: Array{Float64, 2},
    A :: Array{Int, 2},
    b :: Vector{Int}
)
    order.cost_matrix = build_order(new_order, A, b)
end

end
