module GBTools

using LinearAlgebra: I

function has_slacks(A :: Matrix{Int}) :: Bool
    m, n = size(A)
    found_slack = fill(false, m)
    slack_indices = Int[]
    for i in 1:m
        for j in (m+1):n
            if A[i, j] == 1
                found_slack[i] = true
                push!(slack_indices, j)
                break
            end
        end
    end
    return all(found_slack) && (length(slack_indices) == length(unique(slack_indices)))
end

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

function diff(
    gb1 :: Vector{Vector{Int}}, 
    gb2 :: Vector{Vector{Int}}
) :: Vector{Vector{Int}}
    missing_elements = Vector{Int}[]
    for g in gb1
        if !(g in gb2)
            push!(missing_elements, g)
        end
    end
    return missing_elements
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

end
