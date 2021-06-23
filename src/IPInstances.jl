"""
TODO implement 'projections'
- 4ti2 implements projections where some variables are set to unrestricted
- I should do this as well, as soon as I finish implementing the above point on
unrestricted variables
"""
module IPInstances

export IPInstance, original_matrix, original_rhs, original_upper_bounds,
    original_objective, nonnegative_vars, invert_permutation, is_nonnegative,
    is_bounded, unboundedness_proof, update_objective!, nonnegativity_relaxation,
    random_instance

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
        return A, b, C, u, nonnegative
    end
    m, n = size(A)
    In = Matrix{Int}(I, n, n)
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
    new_u = [u; [typemax(Int) for _ in 1:(n+m)]]
    new_nonnegative = [nonnegative; [true for _ in 1:(n+m)]] #slacks are non-negative
    return new_A, new_b, new_C, new_u, new_nonnegative
end

"""
Represents an instance of a problem

min C * x
s.t. A * x = b
0 <= x <= u
x in ZZ^n

The instance is stored in normalized form, with permuted variables so that
the variables appear in the following order: bounded, non-negative, unrestricted.
"""
struct IPInstance
    #Problem data
    A :: Array{Int, 2}
    b :: Vector{Int}
    C :: Array{Float64, 2}
    u :: Vector{Int}

    #Data relative to permutation and variable types
    bounded_end :: Int #index of last bounded variable
    nonnegative_end :: Int #index of last non-negative variable
    permutation :: Vector{Int}
    inverse_permutation :: Vector{Int}

    #Problem metadata
    orig_cons :: Int #constraints before normalization
    orig_vars :: Int #variables before normalization
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
        C :: Array{T, 2},
        u :: Vector{Int},
        nonnegative :: Union{Nothing, Vector{Bool}} = nothing;
        apply_normalization :: Bool = true,
        invert_objective :: Bool = true
    ) where {T <: Real}
        m, n = size(A)
        @assert m == length(b)
        @assert n == size(C, 2)
        @assert n == length(u)
        @assert isnothing(nonnegative) || n == length(nonnegative)
        #If no non-negativity constraints are specified, assume all variables
        #are non-negative
        if isnothing(nonnegative)
            nonnegative = [ true for _ in 1:n ]
        end
        #Normalization of the data to the form Ax = b, minimization...
        A, b, C, u, nonnegative = normalize(
            A, b, C, u, nonnegative,
            apply_normalization=apply_normalization,
            invert_objective=invert_objective
        )
        new_m, new_n = size(A)
        C = Float64.(C)
        #Setting up linear relaxation model
        model, model_vars, model_cons = SolverTools.relaxation_model(A, b, C, u, nonnegative)
        #Checks feasibility of the linear relaxation
        @assert SolverTools.is_feasible(model)
        #Compute boundedness of variables using the model
        bounded = bounded_variables(model, model_vars)
        #Compute a permutation of variables of the given instance such that
        #vars appear in order: bounded, non-negative, unrestricted
        permutation, bounded_end, nonnegative_end = compute_permutation(bounded, nonnegative)
        inverse_perm = invperm(permutation)
        #Permute columns of problem data
        A = A[:, permutation]
        C = C[:, permutation]
        u = u[permutation]
        #Create the normalized instance
        new(A, b, C, u,
            bounded_end, nonnegative_end, permutation, inverse_perm,
            m, n, new_m, new_n, true,
            model, model_vars, model_cons
        )
    end
end

"""
Returns a new IPInstance corresponding to the relaxation of non-negative
constraints given by `nonnegative` in `instance`.
"""
function nonnegativity_relaxation(
    instance :: IPInstance,
    nonnegative :: Vector{Bool}
) :: IPInstance
    return IPInstance(
        instance.A, instance.b, instance.C, instance.u, nonnegative,
        apply_normalization=false,
        invert_objective=false
    )
end

"""
Returns true iff all data in instance.A and instance.b is non-negative and
all variables are non-negative.
"""
function is_nonnegative(
    instance :: IPInstance
) :: Bool
    vars_nonneg = instance.nonnegative_end == instance.n
    a_nonneg = all(ai >= 0 for ai in instance.A)
    b_nonneg = all(bi >= 0 for bi in instance.b)
    return vars_nonneg && a_nonneg && b_nonneg
end

"""
Updates the objective function of `instance` to maximizing its i-th
variable (= minimizing -x_i)
"""
function update_objective!(
    instance :: IPInstance,
    i :: Int
)
    for j in 1:instance.n
        if j == i
            instance.C[1, j] = -1
        else
            instance.C[1, j] = 0
        end
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
# Functions to deal with bounded and non-negative variables
#

"""
Computes a boolean array indicating whether a given variable is bounded.

A variable x_i is bounded for the given model iff
max {x_i | x feasible for model}
is bounded. The model is assumed to be feasible.
"""
function bounded_variables(
    model :: JuMP.Model,
    model_vars :: Vector{JuMP.VariableRef}
) :: Vector{Bool}
    n = length(model_vars)
    bounded = zeros(Bool, n)
    for i in 1:n
        bounded[i] = SolverTools.is_bounded(i, model, model_vars)
    end
    return bounded
end

"""
Returns true iff variable i is nonnegative.
"""
function is_nonnegative(
    i :: Int,
    instance :: IPInstance
) :: Bool
    return i <= instance.nonnegative_end
end

"""
Returns true iff variable i is bounded.
"""
function is_bounded(
    i :: Int,
    instance :: IPInstance
) :: Bool
    return i <= instance.bounded_end
end

"""
Returns a vector indicating whether each variable of `instance` is non-negative.
"""
function nonnegative_vars(
    instance :: IPInstance
) :: Vector{Bool}
    return [ is_nonnegative(i, instance) for i in 1:instance.n ]
end

"""
Returns a vector u in kernel(instance.A) proving that x_i is unbounded.
"""
function unboundedness_proof(
    instance :: IPInstance,
    nonnegative :: Vector{Bool},
    i :: Int
) :: Vector{Int}
    @assert !is_bounded(i, instance)
    model, vars, _ = SolverTools.unboundedness_ip_model(instance.A, nonnegative, i)
    JuMP.optimize!(model)
    if JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        error("Unboundedness model should be feasible.")
    end
    u = Int.(round.(JuMP.value.(vars)))
    return u
end

#
# Permutation-related functions
#

"""
Computes a permutation of variables that puts the variables in the following
order: bounded, then unbounded (but restricted) and finally unrestricted.
This operation should be stable with respect to the initial ordering of
variables. In addition to the permutation, returns the index of the last
bounded variable and the index of the last non-negative variable.

A permutation is represented by a vector perm such that perm[i] = j means
that variable i is sent to j by the permutation.
"""
function compute_permutation(
    bounded :: Vector{Bool},
    nonnegative :: Vector{Bool}
) :: Tuple{Vector{Int}, Int, Int}
    #This code is repetitive and inefficient, but I think it is clear
    #For clarity, I'll leave it this way. It is unlikely to become a bottleneck
    @assert length(bounded) == length(nonnegative)
    n = length(bounded)
    permutation = zeros(Int, n)
    #Set bounded variables first after permutation
    n_bounded = 0
    for i in 1:n
        if bounded[i]
            n_bounded += 1
            permutation[i] = n_bounded
        end
    end
    #Next, we need to set non-negative, unbounded variables
    n_nonnegative = 0
    n_unrestricted = 0
    if n_bounded < n #If every variable is bounded, skip this
        for i in 1:n
            if !bounded[i] && nonnegative[i]
                n_nonnegative += 1
                permutation[i] = n_bounded + n_nonnegative
            end
        end
    end
    #Set unrestricted variables last in the permutation
    if n_bounded + n_nonnegative < n
        for i in 1:n
            if !nonnegative[i]
                n_unrestricted += 1
                permutation[i] = n_bounded + n_nonnegative + n_unrestricted
            end
        end
    end
    @assert n_bounded + n_nonnegative + n_unrestricted == n
    bounded_end = n_bounded
    nonnegative_end = n_bounded + n_nonnegative
    return permutation, bounded_end, nonnegative_end
end

"""
Applies `permutation` to each vector in `vector_set`.
"""
function apply_permutation(
    vector_set :: Vector{Vector{Int}},
    permutation :: Vector{Int}
) :: Vector{Vector{Int}}
    permuted_set = Vector{Int}[]
    for v in vector_set
        @assert length(v) == length(permutation)
        permuted_v = v[permutation] #Julia operation for applying permutation
        push!(permuted_set, permuted_v)
    end
    return permuted_set
end

"""
Inverts the variable permutation of `instance` over `vector_set`.
"""
function invert_permutation(
    vector_set :: Vector{Vector{Int}},
    instance :: IPInstance
) :: Vector{Vector{Int}}
    return apply_permutation(vector_set, instance.inverse_permutation)
end

#
# Generating some IPInstances
#

"""
Returns a random feasible IPInstance with `m` constraints and `n` variables.
"""
function random_instance(
    m :: Int,
    n :: Int
) :: IPInstance
    instance = nothing
    feasible = false
    while !feasible
        #Build random instance in these parameters
        A = rand(-5:5, m, n)
        b = rand(5:20, m)
        C = rand(-10:-1, 1, n)
        u = [ typemax(Int) for _ in 1:n ]
        instance = IPInstance(A, b, C, u)
        #Check feasibility
        model, _, _ = SolverTools.feasibility_model(
            instance.A, instance.b, instance.u, nonnegative_vars(instance), Int
        )
        feasible = SolverTools.is_feasible(model)
    end
    return instance
end

end
