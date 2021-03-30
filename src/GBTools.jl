module GBTools

using LinearAlgebra: I
using JuMP
using Clp

function isincluded(
    gb1 :: Vector{Vector{Int}},
    gb2 :: Vector{Vector{Int}}
) :: Bool
    for g in gb1
        if !(g in gb2)
            return false
        end
    end
    return true
end

function isequal(
    gb1 :: Vector{Vector{Int}},
    gb2 :: Vector{Vector{Int}}
) :: Bool
    return isincluded(gb1, gb2) && isincluded(gb2, gb1)
end

function tomatrix(
    gb :: Vector{Vector{Int}}
) :: Array{Int, 2}
    M = foldl(hcat, gb)
    return M'
end

function tovector(
    gb :: Array{Int, 2}
) :: Vector{Vector{Int}}
    return [ gb[i, :] for i in 1:size(gb, 1) ]
end

"""
Transforms a problem in the form:
max C * x
s.t. Ax <= b
0 <= x <= u

to something of the form
max C * x
s.t. Ax == b
x == u

by adding slack variables.
"""
function normalize(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{T, 2},
    u :: Vector{Int};
    apply_normalization :: Bool = true,
    invert_objective :: Bool = true
) :: Tuple{Array{Int, 2}, Vector{Int}, Array{Float64, 2}, Vector{Int}} where {T <: Real}
    if !apply_normalization
        return A, b, C, u
    end
    m, n = size(A)
    In = Matrix{Int}(I, n, n)
    Zn = zeros(Int, n, n)
    Im = Matrix{Int}(I, m, m)
    Znm = zeros(Int, n, m)
    Zmn = zeros(Int, m, n)
    new_A = [A Im Zmn; In Znm In]
    new_b = [b; u]
    #The reductions without fullfilter only work correctly if the problem
    #is in minimization form. Thus we take the opposite of C instead, as
    #this is easier than changing everything else
    sign = invert_objective ? -1 : 1
    new_C = [sign * C zeros(Int, size(C, 1), n + m)]
    new_u = [u; [typemax(Int) for i in 1:(n+m)]]
    return new_A, new_b, new_C, new_u
end

"""
Returns a matrix representing the grevlex order for `n` variables with
x_n > x_{n-1} > ... > x_1
"""
function grevlex_matrix(
    n :: Int
) :: Array{Int, 2}
    grevlex = Array{Int, 2}(undef, n, n)
    for i in 1:n
        for j in 1:n
            grevlex[i, j] = i <= j ? 1 : 0
        end
    end
    return grevlex
end

"""
Returns a matrix representing the lex order for `n` variables with
x_1 > x_2 > ... > x_n.
"""
function lex_matrix(
    n :: Int
) :: Array{Int, 2}
    return Matrix{Int}(I, n, n)
end

"""
Returns a matrix representing the lex order for `n` variables with
x_1 < x_2 < ... < x_n.

This is the tiebreaking order used by 4ti2.
"""
function revlex_matrix(
    n :: Int
) :: Array{Int, 2}
    return Matrix{Int}(-I, n, n)
end

"""
Computes a strictly positive vector in the row span of `A` using
linear programming. Assumes Ax = b is feasible, and
max x
s.t. Ax = b
x >= 0
is bounded. Given these conditions, the dual variables of the
constraints of the above LP give a positive row span vector.

TODO Actually this thing is giving me some zeroes... check
"""
function positive_row_span(
    A :: Array{Int, 2},
    b :: Vector{Int}
) :: Vector{Float64}
    m, n = size(A)
    model = Model(Clp.Optimizer)
    set_silent(model)
    @variable(model, x[1:n] >= 0)
    constraints = []
    for i in 1:m
        ai = A[i, :]
        con = @constraint(model, ai' * x == b[i])
        push!(constraints, con)
    end
    @objective(model, Max, sum(x[i] for i in 1:n))
    optimize!(model)
    #TODO Check if duals are available, etc
    return A' * shadow_price.(constraints)
end

end
