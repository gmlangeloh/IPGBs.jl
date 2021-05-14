"""
- TODO compute and store which variables are bounded and which aren't
- TODO store the unrestricted variables (useful for Markov basis stuff)
"""
module IPInstances

export IPInstance, original_matrix, original_rhs, original_upper_bounds,
    original_objective

import LinearAlgebra: I

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
Represents an instance of a problem

min C * x
s.t. A * x = b
0 <= x <= u
x in ZZ^n
"""
struct IPInstance
    A :: Array{Int, 2}
    b :: Vector{Int}
    C :: Array{Float64, 2}
    u :: Vector{Int}
    orig_cons :: Int
    orig_vars :: Int
    m :: Int
    n :: Int
    sense :: Bool #true if minimization

    #TODO put a parameter to determine whether it is minimization or not
    function IPInstance(A, b, C, u)
        m, n = size(A)
        @assert m == length(b)
        @assert n == size(C, 2)
        @assert n == length(u)
        A, b, C, u = normalize(
            A, b, C, u, apply_normalization=true, invert_objective=true
        )
        new_m, new_n = size(A)
        C = Float64.(C)
        new(A, b, C, u, m, n, new_m, new_n, true)
    end
end

function original_matrix(
    instance :: IPInstance
) :: Array{Int, 2}
    m = instance.orig_cons
    n = instance.orig_vars
    return instance.A[1:m, 1:n]
end

function original_rhs(
    instance :: IPInstance
) :: Vector{Int}
    return instance.b[1:instance.orig_cons]
end

function original_upper_bounds(
    instance :: IPInstance
) :: Vector{Int}
    return instance.u[1:instance.orig_vars]
end

function original_objective(
    instance :: IPInstance
) :: Array{Float64, 2}
    return instance.C[:, 1:instance.orig_vars]
end

end
