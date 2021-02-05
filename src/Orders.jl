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
    C :: Array{Int, 2}
) :: Array{Int, 2}
    m = size(C, 1)
    n = size(C, 2)
    if m < n #Insufficient tiebreaking in C
        tiebreaker = GBTools.grevlex_matrix(n)
        full_matrix = vcat(C, tiebreaker)
    else
        full_matrix = C
    end
    #Store the transpose to exploit the fact Julia stores matrices
    #in column-major form. This helps performance of cmp and lt.
    return full_matrix'
end

"""
Implements a monomial order from a given cost matrix, including a grevlex
tiebreaker if necessary.

- TODO support other tiebreakers such as lex, as well as other permutations
of the grevlex variables
"""
mutable struct MonomialOrder <: GBOrder
    cost_matrix :: Array{Int, 2}
    #Probably should carry a tiebreaker around too
    MonomialOrder(costs :: Array{Int, 2}) = new(build_order(costs))
end

#TODO should the length be the number of rows or columns?
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
function is_inverted(
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
Changes the order stored in MonomialOrder to a new order. Rebuilds
tiebreaker if necessary.
"""
function change_ordering!(
    order :: MonomialOrder,
    new_order :: Array{Int, 2}
)
    order.cost_matrix = build_order(new_order)
end

end
