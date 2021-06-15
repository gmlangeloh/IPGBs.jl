"""
TODO implement 'projections'
- 4ti2 implements projections where some variables are set to unrestricted
- I should do this as well, as soon as I finish implementing the above point on
unrestricted variables
"""
module IPInstances

export IPInstance, original_matrix, original_rhs, original_upper_bounds,
    original_objective

import LinearAlgebra: I
import JuMP

using IPGBs.SolverTools

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
    u :: Vector{Int},
    nonnegative :: Vector{Bool};
    apply_normalization :: Bool = true,
    invert_objective :: Bool = true
) :: Tuple{Array{Int, 2}, Vector{Int}, Array{Float64, 2}, Vector{Int}, Vector{Bool}} where {T <: Real}
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
    new_nonnegative = [nonnegative; [true for i in 1:(n+m)]] #slacks are non-negative
    return new_A, new_b, new_C, new_u, new_nonnegative
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
    nonnegative :: Vector{Bool} #true iff the i-th var has a non-negativity constr
    orig_cons :: Int
    orig_vars :: Int
    m :: Int
    n :: Int
    sense :: Bool #true if minimization

    #Store a linear relaxation of this instance as a JuMP model
    #It is used to check whether variables are bounded
    model :: JuMP.Model
    model_vars :: Vector{JuMP.VariableRef}
    model_cons :: Vector{JuMP.ConstraintRef} #TODO not a concrete type, fix this

    #TODO put a parameter to determine whether it is minimization or not
    function IPInstance(
        A :: Array{Int, 2},
        b :: Vector{Int},
        C :: Array{Int, 2},
        u :: Vector{Int},
        nonnegative :: Union{Nothing, Vector{Bool}} = nothing
    )
        m, n = size(A)
        @assert m == length(b)
        @assert n == size(C, 2)
        @assert n == length(u)
        @assert isnothing(nonnegative) || n == length(nonnegative)
        #If no non-negativity constraints are specified, assume all variables
        #are non-negative
        if isnothing(nonnegative)
            nonnegative = [ true for i in 1:n ]
        end
        #Normalization of the data to the form Ax = b, minimization...
        A, b, C, u, nonnegative = normalize(
            A, b, C, u, nonnegative, apply_normalization=true, invert_objective=true
        )
        new_m, new_n = size(A)
        C = Float64.(C)
        #Setting up linear relaxation model
        model, model_vars, model_cons = SolverTools.relaxation_model(A, b, C, u, nonnegative)
        #Checks feasibility of the linear relaxation
        @assert SolverTools.is_feasible(model)
        new(A, b, C, u, nonnegative, m, n, new_m, new_n, true, model, model_vars, model_cons)
    end
end

#
# The following functions are used to obtain original, non-normalized data
#

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

#
# Functions to deal with boundedness of components / variables
# A variable is x_i is bounded iff max {x_i \mid x feasible} is bounded.
# This is equivalent to the linear relaxation of this problem being bounded
#

function is_bounded(
    i :: Int,
    instance :: IPInstance
) :: Bool
    return SolverTools.is_bounded(i, instance.model, instance.model_vars)
end

end
