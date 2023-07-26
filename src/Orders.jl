module Orders

export GBOrder, MonomialOrder, is_inverted, order_costs

import LinearAlgebra: I

using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances
using MIPMatrixTools.SolverTools

using IPGBs

"""
This is specialized by MonomialOrder in case of Buchberger's algorithm and by
ModuleMonomialOrder in case of Signature-based algorithms.
"""
abstract type GBOrder <: Base.Ordering end

is_inverted(:: GBOrder, :: AbstractVector{Int}) = error("Not implemented.")
order_costs(:: GBOrder, :: AbstractVector{Int}) = error("Not implemented.")

"""
    tiebreak(C :: Matrix{Float64}) :: Matrix{Float64}

Return a matrix corresponding to the same monomial order as `C`, tiebreaking
it using the revlex order if `C` does not already specify a full monomial order
itself.
"""
function tiebreak(C :: Matrix{Float64})
    m, n = size(C)
    if m < n #Insufficient tiebreaking in C
        tiebreaker = GBTools.revlex_matrix(n)
        full_matrix = vcat(C, tiebreaker)
    else
        full_matrix = C
    end
    return full_matrix
end

"""
    positive_first_row!(C :: Matrix{Float64}, A :: Matrix{Int}, b :: Vector{Int})

Modify `C` so that its first row is strictly positive using dual data from
the LP min C*x s.t. Ax = b, x >= 0.

Assumes the given LP is bounded and feasible.
"""
function positive_first_row!(
    C :: Matrix{Float64},
    A :: Matrix{Int},
    b :: Vector{Int}
)
    n = size(C, 2)
    @assert n == size(A, 2)
    if any(C[1, j] < 0 for j in 1:n)
        d = SolverTools.positive_row_span(A, b)
        if !isnothing(d)
            lambda = 1 + maximum(-C[1, j] / d[j] for j in 1:n)
            for j in 1:n
                C[1, j] += lambda * d[j]
            end
        end
    end
end

"""
    projection_order(C :: Matrix{Float64}, A :: Matrix{Int}, b :: Vector{Int}, num_vars :: Int)

Return a matrix corresponding to the same monomial order as `C` but using only
the first `num_vars` variables.

Uses the dual of the LP min C*x s.t. Ax = b, x >= 0 to find such an order.
"""
function projection_order(
    C :: Matrix{Float64},
    A :: Matrix{Int},
    b :: Vector{Int},
    num_vars :: Int
) :: Matrix{Float64}
    s = SolverTools.optimal_row_span(A, b, C)
    projected_obj = copy(C)
    projected_obj[1, :] -= s
    #TODO This assertion may fail when the LP above is degenerate
    #I should prove that I can project away the objective regardless
    #@assert(all(IPGBs.is_approx_zero(projected_obj[1, i]) for i in (num_vars+1):size(C, 2)))
    return projected_obj[:, 1:num_vars]
end

"""
    normalize_order(C :: Matrix{Float64}, A :: Matrix{Int}, b :: Vector{Int}, num_vars :: Int)

Return a cost matrix in column-major form giving the same monomial order as `C`
using only the first num_vars variables.

Normalization uses the dual data of the LP min C*x s.t. Ax = b, x >= 0.
"""
function normalize_order(
    C::Matrix{Float64},
    A::Matrix{Int},
    b::Vector{Int},
    num_vars::Int
)::Matrix{Float64}
    cost_matrix = zeros(Float64, num_vars, num_vars)
    if num_vars == size(A, 2)
        cost_matrix = tiebreak(C)
        positive_first_row!(cost_matrix, A, b)
    else
        #Compute order projected into the first num_vars
        proj = projection_order(C, A, b, num_vars)
        cost_matrix = tiebreak(proj)
    end
    #Store the transpose to exploit the fact Julia stores matrices
    #in column-major form. This helps performance of cmp and lt.
    return cost_matrix'
end

"""
Implements a monomial order from a given cost matrix, including a grevlex
tiebreaker if necessary.

The cost matrix is represented in column-major order to speed up comparisons.
"""
mutable struct MonomialOrder <: GBOrder
    cost_matrix :: Array{Float64, 2}
    tiebreaker :: Symbol
    is_minimization :: Bool
    num_costs :: Int
end

function MonomialOrder(
    costs::Matrix{Float64},
    A::Matrix{Int},
    b::Vector{Int},
    is_minimization::Bool,
    num_vars :: Union{Nothing, Int} = nothing
)
    num_costs = size(costs, 1)
    max_vars = size(A, 2)
    if !isnothing(num_vars)
        max_vars = num_vars
    end
    return MonomialOrder(normalize_order(costs, A, b, max_vars), :invlex,
                         is_minimization, num_costs)
end

MonomialOrder(C :: Matrix{S}, A, b, min, nvars = nothing) where {S <: Real} =
    MonomialOrder(Float64.(C), A, b, min, nvars)

function MonomialOrder(
    instance :: IPInstance,
    num_vars :: Union{Nothing, Int} = nothing
)
    return MonomialOrder(instance.C, instance.A, instance.b, true, num_vars)
end

"""
    lex_order(n :: Int) :: MonomialOrder

Return the lex monomial order with `x1 > x2 > ... > xn`.
"""
function lex_order(
    n :: Int
) :: MonomialOrder
    C = Matrix{Float64}(I(n))
    #No normalization is necessary for a lex order.
    #Note: column-major or row-major is irrelevant here, as transpose(I) == I.
    return MonomialOrder(C, :lex, true)
end

"""
    order_costs(
    order :: MonomialOrder,
    v :: AbstractVector{Int}
) :: Vector{Int}

The cost of some vector according to the cost matrix of `order`.

This is the vector weighted by the first column of the cost matrix.
"""
function order_costs(
    order :: MonomialOrder,
    v :: AbstractVector{Int}
) :: Vector{Int}
    return [ round(Int, order.cost_matrix[:, j]' * v) 
        for j in 1:order.num_costs]
end

Base.length(o :: MonomialOrder) = size(o.cost_matrix, 2)

#This is implemented this way for locality and performance. It is only used
#in cmp / lt below regardless.
Base.getindex(o :: MonomialOrder, i, j) = o.cost_matrix[j, i]

function Base.cmp(
    order :: MonomialOrder,
    v1 :: T,
    v2 :: T
) :: Int where {T <: AbstractVector{Int}}
    @assert length(v1) == length(v2)
    #TODO check if I should really change the sign in the maximization case
    minimization = order.is_minimization ? 1 : 1
    for i in 1:length(order)
        #Compute o' * (v1 - v2) without allocations
        s = 0.0
        for j in 1:length(v1)
            s += order[i, j] * (v1[j] - v2[j])
        end
        if isapprox(s, 0.0, atol=1e-8)
            continue
        elseif s < 0.0
            return -1 * minimization
        elseif s > 0.0
            return 1 * minimization
        end
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
    invert_maximization(order :: MonomialOrder, result :: Bool)

Invert `result` when `order` is a maximization order.
"""
function invert_maximization(
    _ :: MonomialOrder,
    result :: Bool
) :: Bool
    #TODO I don't know if this is correct in the maximization case...
    #if order.is_minimization
    #    return result
    #end
    #return !result
    return result
end

"""
    is_inverted_generic(order :: MonomialOrder, v :: T) :: Bool where {T <: AbstractVector{Int}}

Return true iff the trailing and leading terms of `v` are inverted with
respect to `order`.
"""
function is_inverted_generic(
    order :: MonomialOrder,
    v :: T
) :: Bool where {T <: AbstractVector{Int}}
    for i in 1:length(order)
        s = 0.0
        for j in 1:length(v)
            s += order[i, j] * v[j]
        end
        if isapprox(s, 0.0, atol=1e-8)
            #Avoids numeric issues
            continue
        elseif s < 0.0
            return invert_maximization(order, true)
        elseif s > 0.0
            return invert_maximization(order, false)
        end
    end
    return false #This really should not happen
end

"""
    is_inverted_invlex(order :: MonomialOrder, v :: T) :: Bool where {T <: AbstractVector{Int}}

Return true iff the trailing and leading terms of `v` are inverted with
respect to the invlex order.

`v` is inverted according to invlex iff its first non-zero entry is positive.
Inversion testing is thus done more efficiently in the invlex case.
"""
function is_inverted_invlex(
    order :: MonomialOrder,
    v :: T
) :: Bool where {T <: AbstractVector{Int}}
    i = 1
    while i <= length(v) && v[i] == 0
        i += 1
    end
    if i > length(v)
        return invert_maximization(order, false)
    end
    return invert_maximization(order, v[i] > 0)
end

"""
    is_inverted(
    order :: MonomialOrder,
    v :: T,
    costs :: Vector{Int}
) :: Bool where {T <: AbstractVector{Int}}

Efficiently return whether v's leading and trailing terms are inverted.
"""
function is_inverted(
    order :: MonomialOrder,
    v :: T,
    costs :: AbstractVector{Int}
) :: Bool where {T <: AbstractVector{Int}}
    for c in costs
        if c < 0
            return invert_maximization(order, true)
        elseif c > 0
            return invert_maximization(order, false)
        end
    end
    if order.tiebreaker == :invlex
        return is_inverted_invlex(order, v)
    end
    return is_inverted_generic(order, v)
end

"""
    change_ordering!(order :: MonomialOrder, new_order :: Matrix{Float64}, A :: Matrix{Int}, b :: Vector{Int})

Change the order stored in MonomialOrder to a new order.

The tiebreaker is rebuilt if necessary.
"""
function change_ordering!(
    order :: MonomialOrder,
    new_order :: Array{Float64, 2},
    A :: Array{Int, 2},
    b :: Vector{Int}
)
    order.cost_matrix = normalize_order(new_order, A, b, size(A, 2))
end

end
