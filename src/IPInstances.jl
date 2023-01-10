module IPInstances

export IPInstance, nonnegative_vars, is_bounded, unboundedness_proof, update_objective!, nonnegativity_relaxation, group_relaxation, lift_vector, truncation_weight

import LinearAlgebra: I
using AbstractAlgebra
using JuMP

using IPGBs
using IPGBs.SolverTools

const AlgebraInt = AbstractAlgebra.Integers{Int}()

"""
    normalize_ip(A :: Matrix{Int}, b :: Vector{Int}, c :: Matrix{T}, u :: Vector{<: Union{Int, Nothing}}, nonnegative :: Vector{Bool}; ...) where {T <: Real}

Transform a problem in the form:
max C * x
s.t. Ax <= b
0 <= x <= u

to something of the form
max C * x
s.t. Ax == b
x == u

by adding slack variables.
"""
function normalize_ip(
    A::Matrix{Int},
    b::Vector{Int},
    C::Matrix{T},
    u::Vector{<: Union{Int,Nothing}},
    nonnegative::Vector{Bool};
    apply_normalization::Bool = true,
    invert_objective::Bool = true
)::Tuple{Array{Int,2},Vector{Int},Array{Float64,2},Vector{Union{Int,Nothing}},Vector{Bool}} where {T<:Real}
    if !apply_normalization
        return A, b, C, u, nonnegative
    end
    m, n = size(A)
    k = count(!isnothing(u[i]) for i in 1:length(u))
    Ik = Matrix{Int}(I, k, k)
    UBkn = zeros(Int, k, n)
    ubs = Int[]
    i = 1
    for j in 1:n
        if !isnothing(u[i])
            push!(ubs, u[i])
            UBkn[i, j] = 1
            i += 1
        end
    end
    Im = Matrix{Int}(I, m, m)
    Zkm = zeros(Int, k, m)
    Zmk = zeros(Int, m, k)
    new_A = [A Im Zmk; UBkn Zkm Ik]
    new_b = [b; ubs]
    #The reductions without fullfilter only work correctly if the problem
    #is in minimization form. Thus we take the opposite of C instead, as
    #this is easier than changing everything else
    sign = invert_objective ? -1 : 1
    new_C = [sign * C zeros(Int, size(C, 1), k + m)]
    new_u = Union{Int, Nothing}[]
    for val in u
        push!(new_u, val)
    end
    for _ in 1:(k+m)
        push!(new_u, nothing)
    end
    new_nonnegative = [nonnegative; [true for _ in 1:(k+m)]] #slacks are non-negative
    return new_A, new_b, new_C, new_u, new_nonnegative
end

"""
    hnf_lattice_basis(A :: Matrix{Int})

    Return a row basis for the lattice ker(A) computed using the Upper
    Hermite Normal Form. The entries of this basis tend to be smaller than
    those computed directly from kernel(A).
"""
function hnf_lattice_basis(A :: Matrix{Int})
    m, n = size(A)
    mat_A = matrix(AlgebraInt, transpose(A))
    r = rank(mat_A)
    #Transpose and append identity matrix, so that the lattice basis appears
    #as the last few rows / columns of the uhnf.
    tA = hcat(mat_A, identity_matrix(AlgebraInt, n))
    #tA is a n x (m + n) matrix.
    #hnf_cohen is often slightly faster than hnf
    I = identity_matrix(tA, n)
    #Even though there are apparently no guarantees, running hnf_kb! over 64-bit
    #ints does work. Running hnf_cohen! here instead doesn't, though.
    AbstractAlgebra.hnf_kb!(tA, I)
    #The basis is in the last few rows and columns of H
    basis = tA[(r+1):n, (m+1):(n+m)]
    return basis, r #Row basis of the lattice
end

"""
    fiber_solution(A :: Matrix{Int}, b :: Vector{Int}) :: Vector{Int}

    A solution to Ax = b, that is, an element of the fiber of right-hand side
    `b` in the lattice ker(A).
"""
function fiber_solution(A :: Matrix{Int}, b :: Vector{Int}) :: Vector{Int}
    m, n = size(A)
    mat_A = matrix(AlgebraInt, A)
    mat_b = matrix(AlgebraInt, m, 1, b)
    x = solve(mat_A, mat_b)
    return Int.(reshape(Array(x), n))
end

"""
Represents an instance of a problem

min C * x

s.t. A * x = b

0 <= x <= u

x in ZZ^n

The instance is stored in normalized form, with permuted variables so that
the variables appear in the following order: bounded, non-negative but unbounded, unrestricted.
"""
struct IPInstance
    #Problem data
    A :: Array{Int, 2}
    b :: Vector{Int}
    C :: Array{Float64, 2}
    u :: Vector{Union{Int, Nothing}}

    #Data relative to permutation and variable types
    bounded_end :: Int #index of last bounded variable
    nonnegative_end :: Int #index of last non-negative variable
    permutation :: Vector{Int}
    inverse_permutation :: Vector{Int}

    #Problem metadata
    orig_cons :: Int #constraints before normalization
    orig_vars :: Int #variables before normalization
    m :: Int #number of constraints after normalization
    n :: Int #number of variables after normalization
    sense :: Bool #true if minimization

    #Store a linear relaxation of this instance as a JuMP model
    #It is used to check whether variables are bounded
    model :: JuMP.Model
    model_vars :: Vector{JuMP.VariableRef}
    model_cons :: Vector{JuMP.ConstraintRef} #TODO: not a concrete type, fix this

    #Lattice-related information
    lattice_basis :: Generic.MatSpaceElem{Int} #Row basis
    rank :: Int
    fiber_solution :: Vector{Int} #v such that Av = b. Not necessarily non-negative.
    originally_bounded :: Vector{Bool}

    #TODO: put a parameter to determine whether it is minimization or not
    function IPInstance(
        A::Array{Int,2},
        b::Vector{Int},
        C::Array{T,2},
        u::Vector{<:Union{Int,Nothing}},
        nonnegative::Union{Nothing,Vector{Bool}} = nothing;
        apply_normalization::Bool = true,
        invert_objective::Bool = true
    ) where {T<:Real}
        m, n = size(A)
        @assert m == length(b)
        @assert n == size(C, 2)
        @assert n == length(u)
        @assert isnothing(nonnegative) || n == length(nonnegative)
        #If no non-negativity constraints are specified, assume all variables
        #are non-negative
        if isnothing(nonnegative)
            nonnegative = [true for _ in 1:n]
        end
        #Normalization of the data to the form Ax = b, minimization...
        A, b, C, u, nonnegative = normalize_ip(
            A, b, C, u, nonnegative,
            apply_normalization = apply_normalization,
            invert_objective = invert_objective
        )
        new_m, new_n = size(A)
        C = Float64.(C)
        #Create a JuMP model to compute bounded variables
        model, model_vars, model_cons = SolverTools.relaxation_model(A, b, C, u, nonnegative)
        #Compute a permutation of variables of the given instance such that
        #vars appear in order: bounded, non-negative, unrestricted
        bounded = bounded_variables(model, model_vars, model_cons, b)
        permutation, bounded_end, nonnegative_end = compute_permutation(bounded, nonnegative)
        inverse_perm = invperm(permutation)
        #Permute columns of problem data
        A = A[:, permutation]
        C = C[:, permutation]
        u = u[permutation]
        #Update the JuMP model with the permutation info
        model, model_vars, model_cons = SolverTools.relaxation_model(A, b, C, u, nonnegative)
        #Checks feasibility of the linear relaxation
        @assert SolverTools.is_feasible(model)
        #Compute boundedness of variables using the model
        SolverTools.set_jump_objective!(model, :Min, C[1, :])
        #Compute lattice information
        basis, rnk = hnf_lattice_basis(A)
        fiber_sol = fiber_solution(A, b)
        #Create the normalized instance
        new(A, b, C, u,
            bounded_end, nonnegative_end, permutation, inverse_perm,
            m, n, new_m, new_n, true,
            model, model_vars, model_cons,
            basis, rnk, fiber_sol, bounded
        )
    end
end

function Base.show(io::IO, instance::IPInstance)
    obj = "min $(instance.C) \n"
    constr = ""
    for i in 1:size(instance.A, 1)
        row = instance.A[i, :]
        str = "$row = $(instance.b[i]) \n"
        constr *= str
    end
    bounds = ""
    for i in 1:length(instance.u)
        if isnothing(instance.u[i])
            continue
        end
        bounds *= "0 <= x$(i) <= $(instance.u[i])"
        if i < length(instance.u)
            bounds *= "\n"
        end
    end
    final = obj * constr * bounds
    print(io, final)
end

"""
    extract_constraint(model :: JuMP.Model, c :: JuMP.ConstraintRef, x :: Vector{JuMP.VariableRef})

Extract numerical coefficients from a JuMP constraint, returning a vector
with left-hand side coefficients and the right-hand side value.

Assumes this is a scalar constraint (= a single constraint).
"""
function extract_constraint(
    model :: JuMP.Model,
    c :: JuMP.ConstraintRef,
    x :: Vector{JuMP.VariableRef}
)
    lhs_data = MOI.get(model, MOI.ConstraintFunction(), c)
    #List of pairs (coef, variable) where variable is a MOI.VariableIndex
    #this index can be compared to index(:: VariableRef)
    coef_vars = [(term.coefficient, term.variable) for term in lhs_data.terms]
    #Extract row of the constraint matrix, write to a
    a = zeros(Int, length(x))
    for (coef, var_index) in coef_vars
        for j in 1:length(x)
            if index(x[j]) == var_index
                a[j] = coef
                break
            end
        end
    end
    rhs_data = MOI.get(model, MOI.ConstraintSet(), c)
    if hasfield(typeof(rhs_data), :value)
        b = rhs_data.value
    elseif hasfield(typeof(rhs_data), :lower)
        b = rhs_data.lower
    elseif hasfield(typeof(rhs_data), :upper)
        b = rhs_data.upper
    else
        error("Unknown JuMP RHS type: ", typeof(rhs_data))
    end
    return a, b
end

"""
    extract_bound(model :: JuMP.Model, c :: JuMP.ConstraintRef, x :: Vector{JuMP.VariableRef})

Extract lower / upper bound value from a VariableRef type constraint `c`.
"""
function extract_bound(
    model :: JuMP.Model,
    c :: JuMP.ConstraintRef,
    x :: Vector{JuMP.VariableRef}
)
    #Upper and lower bounds have to be turned into explicit constraints,
    #except for 0 lower bounds.
    var_index = MOI.get(model, MOI.ConstraintFunction(), c)
    x_index = 0
    for j in 1:length(x)
        if index(x[j]) == var_index
            x_index = j
            break
        end
    end
    wrapped_bound = MOI.get(model, MOI.ConstraintSet(), c)
    bound_value = 0
    if typeof(wrapped_bound) <: MOI.GreaterThan #Lower bound
        bound_value = wrapped_bound.lower
    elseif typeof(wrapped_bound) <: MOI.LessThan #Upper bound
        bound_value = wrapped_bound.upper
    else
        error("Unknown variable bound type, " * string(typeof(wrapped_bound)))
    end
    return x_index, bound_value
end

function extract_objective(obj :: JuMP.AffExpr, x :: Vector{JuMP.VariableRef})
    c = zeros(Int, length(x))
    for j in 1:length(x)
        coef = obj.terms[x[j]]
        c[j] = Int(round(coef))
    end
    return c
end

function extract_objective(obj::JuMP.VariableRef, x::Vector{JuMP.VariableRef})
    # JuMP represents the objective function as a single VariableRef when
    #possible.
    # This case needs to be treated separately here.
    c = zeros(Int, length(x))
    index = obj.index.value
    c[index] = 1
    return c
end

function extract_objective(
    model :: JuMP.Model,
    x :: Vector{JuMP.VariableRef}
) :: Vector{Int}
    c = extract_objective(objective_function(model), x)
    if objective_sense(model) == MOI.MAX_SENSE
        c = -c
    end
    #If the objective sense is minimization, no normalization of the
    #objective function coefficients is needed.
    return c
end

function IPInstance(model::JuMP.Model)
    #Extract A, b, c from the model.
    n = num_variables(model)
    x = all_variables(model)
    rows = []
    rhs = []
    ineq_directions = []
    upper_bounds = []
    lower_bounds = []
    #Extract all data from the JuMP model
    for (t1, t2) in list_of_constraint_types(model)
        cs = all_constraints(model, t1, t2)
        for constraint in cs
            if t1 <: AffExpr #Linear constraint
                push!(ineq_directions, t2)
                a, b = extract_constraint(model, constraint, x)
                push!(rows, a)
                push!(rhs, b)
            elseif t1 <: VariableRef #Variable upper / lower bound
                if t2 <: MOI.Integer
                    #Explicit integrality constraints are unnecessary
                    continue
                elseif t2 <: MOI.ZeroOne
                    #Binary constraints for variables added as upper bounds
                    push!(upper_bounds, (constraint.index.value, 1))
                    continue
                end
                bound = extract_bound(model, constraint, x)
                if t2 <: MOI.LessThan #Upper bound
                    push!(upper_bounds, bound)
                elseif t2 <: MOI.GreaterThan #Lower bound
                    push!(lower_bounds, bound)
                end
            end
        end
    end
    #Build matrix representation of the IP min{cx | Ax == b, x >= 0}
    c = extract_objective(model, x)
    #Add upper and lower bounds to A, whenever necessary
    for (var, lb) in lower_bounds
        if !IPGBs.is_approx_zero(lb) #Zero lower bounds may be ignored
            #TODO: Zero lbs are important for project-and-lift, add them later
            #separately from the rest of the data
            new_row = zeros(Int, n)
            new_row[var] = 1
            push!(rows, new_row)
            push!(rhs, lb)
            push!(ineq_directions, MOI.GreaterThan{Float64})
        end
    end
    for (var, ub) in upper_bounds
        new_row = zeros(Int, n)
        new_row[var] = 1
        push!(rows, new_row)
        push!(rhs, ub)
        push!(ineq_directions, MOI.LessThan{Float64})
    end
    A = Int.(foldl(vcat, map(row -> row', rows)))
    #Add slack variables for all inequalities to A
    m = size(A, 1)
    for i in 1:m
        if ineq_directions[i] <: MOI.LessThan
            #Slack in the <= case has positive coefficients
            new_col = zeros(Int, m, 1)
            new_col[i] = 1
            A = hcat(A, new_col)
        elseif ineq_directions[i] <: MOI.GreaterThan
            #Slack in the >= case has negative coefficients
            new_col = zeros(Int, m, 1)
            new_col[i] = -1
            A = hcat(A, new_col)
        end
    end
    #Build right hand side vector
    b = zeros(Int, m)
    for i in 1:m
        try
            b[i] = Int(round(rhs[i]))
        catch e
            if isa(e, InexactError)
                b[i] = typemax(Int)
            else
                throw(e)
            end
        end
    end
    #Build upper bound vector
    u = fill(typemax(Int), size(A, 2))
    for (var, ub) in upper_bounds
        u[var] = ub
    end
    #Extend c to the slack variables
    for _ in 1:(size(A, 2) - length(c))
        push!(c, 0)
    end
    #Build the IPInstance object. The matrices here are already normalized,
    #so no additional normalization is necessary
    return IPInstance(A, b, reshape(c, (1, length(c))), u,
                      apply_normalization=false)
end

function integer_objective(
    instance :: IPInstance
) :: Array{Int}
    k, n = size(instance.C)
    integer_C = zeros(Int, k, n)
    #Find lcm of the denominators of instance.C
    denoms = [ denominator(Rational(c)) for c in instance.C ]
    l = lcm(denoms)
    #Create integer objective function
    for i in 1:k
        for j in 1:n
            integer_C[i, j] = Int(round(l * instance.C[i, j]))
        end
    end
    return integer_C
end

function has_variable_bound_constraints(instance :: IPInstance) :: Bool
    cols_not_slacks = instance.n - instance.m
    rows_not_variable_bounds = instance.m - cols_not_slacks
    if rows_not_variable_bounds < cols_not_slacks
        return false
    end
    #Check if the last few constraints form an identity matrix
    for i in (rows_not_variable_bounds+1):instance.m
        for j in 1:cols_not_slacks
            if i - rows_not_variable_bounds == j
                if instance.A[i, j] == 1
                    return false #If not 1 in the diagonal, not an identity
                end
            elseif instance.A[i, j] != 0
                #If not 0 elsewhere, not an identity
                return false
            end
        end
    end
    return true
end

function nonnegative_variables(instance :: IPInstance) :: Vector{Bool}
    return [ i <= instance.nonnegative_end for i in 1:instance.n ]
end

"""
    nonnegativity_relaxation(instance :: IPInstance, nonnegative :: Vector{Bool}) :: IPInstance

Return a new IPInstance corresponding to the relaxation of `instance`
consisting of only keeping the non-negativity constraints of variables
marked in `nonnegative`.
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
    group_relaxation(instance :: IPInstance) :: IPInstance

Return a new IPInstance corresponding to the relaxation in `instance` of
the non-negativity constraints of the basic variables in the optimal solution
of its linear relaxation.
"""
function group_relaxation(
    instance :: IPInstance
) :: IPInstance
    basis = SolverTools.optimal_basis!(instance.model)
    nonbasics = [ !variable for variable in basis ]
    @assert count(basis) == instance.m
    #Keep only the nonnegativity constraints on the non-basic variables
    return nonnegativity_relaxation(instance, nonbasics)
end

function lattice_basis_projection(
    instance :: IPInstance
)
    # TODO: I need to guarantee the selected columns are linearly
    #independent!!!
    lattice_rank = instance.n - instance.rank
    return instance.lattice_basis[:, 1:lattice_rank]
end

"""
    group_initial_solution(
    instance :: IPInstance
) :: Vector{Int}

Return a solution to instance.A == instance.b, dropping the non-negativity
constraints.
"""
function initial_solution(
    instance :: IPInstance
) :: Vector{Int}
    mat_A = matrix(AlgebraInt, instance.A)
    mat_b = matrix(AlgebraInt, length(instance.b), 1, instance.b)
    x = AbstractAlgebra.solve(mat_A, mat_b)
    return reshape(convert.(Int, Array(x)), size(instance.A, 2))
end

"""
    lift_vector(v :: Vector{Int}, instance :: IPInstance) :: Vector{Int}

Lift `v` from an (extended) group relaxation to the full problem given by
`instance`.
"""
function lift_vector(
    v :: Vector{Int},
    instance :: IPInstance
) :: Vector{Int}
    #Find a lifted vector to full variables
    col_basis = transpose(lattice_basis_projection(instance))
    coefs = AbstractAlgebra.solve(col_basis, matrix(AlgebraInt, length(v), 1, v))
    full_col_basis = transpose(instance.lattice_basis)
    res = full_col_basis * coefs
    return reshape(Array(res), length(res))
end

function truncation_weight(
    instance :: IPInstance
) :: Tuple{Vector{Float64}, Float64}
    A = Array(Int.(instance.lattice_basis))
    b = instance.fiber_solution
    unbounded = map(x -> !x, instance.originally_bounded)
    return SolverTools.optimal_weight_vector(A, b, unbounded)
end

"""
    nonnegative_data_only(instance :: IPInstance) :: Bool

Return true iff all data in `instance.A` and `instance.b` is non-negative and
all variables are non-negative.
"""
function nonnegative_data_only(
    instance :: IPInstance
) :: Bool
    vars_nonneg = instance.nonnegative_end == instance.n
    a_nonneg = all(ai >= 0 for ai in instance.A)
    b_nonneg = all(bi >= 0 for bi in instance.b)
    return vars_nonneg && a_nonneg && b_nonneg
end

"""
    update_objective!(instance :: IPInstance, i :: Int)

Update the objective function of `instance` to maximizing its `i`-th variable
(= minimizing -x_i)
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

function update_objective!(
    instance :: IPInstance,
    i :: Int,
    sigma :: Vector{Int}
)
    c = SolverTools.bounded_objective(instance.A, i, sigma)
    #We take the negative here to normalize the problem to minimization form
    instance.C[1, :] = -c
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
    bounded_variables(model :: JuMP.Model, model_vars :: Vector{JuMP.VariableRef}) :: Vector{Bool}

Return a boolean array indicating whether a given variable is bounded.

A variable x_i is bounded for the given model iff max {x_i | x feasible for model} is bounded for all feasible RHS.
"""
function bounded_variables(
    model :: JuMP.Model,
    model_vars :: Vector{JuMP.VariableRef},
    model_cons :: Vector{JuMP.ConstraintRef},
    b :: Vector{Int}
) :: Vector{Bool}
    n = length(model_vars)
    bounded = zeros(Bool, n)
    #Set the RHS to 0 to check boundedness for all RHS, by duality
    set_normalized_rhs.(model_cons, 0)
    for i in 1:n
        bounded[i] = SolverTools.is_bounded(i, model, model_vars)
    end
    #Set the RHS back to its original value
    set_normalized_rhs.(model_cons, b)
    return bounded
end

"""
    is_nonnegative(i :: Int, instance :: IPInstance) :: Bool

Return true iff the variable of index `i` in `instance` is nonnegative.
"""
function is_nonnegative(
    i :: Int,
    instance :: IPInstance
) :: Bool
    return i <= instance.nonnegative_end
end

"""
    is_bounded(i :: Int, instance :: IPInstance) :: Bool

Return true iff the variable of index `i` in `instance` is bounded.
"""
function is_bounded(
    i :: Int,
    instance :: IPInstance
) :: Bool
    return i <= instance.bounded_end
end

function is_bounded(
    instance :: IPInstance
) :: Bool
    SolverTools.set_jump_objective!(instance.model, :Min, vec(instance.C[1, :]))
    return SolverTools.is_bounded(instance.model)
end

"""
    nonnegative_vars(instance :: IPInstance) :: Vector{Bool}

Return a boolean vector indicating whether each variable of `instance` is non-negative.
"""
function nonnegative_vars(
    instance :: IPInstance
) :: Vector{Bool}
    return [ is_nonnegative(i, instance) for i in 1:instance.n ]
end

"""
    unboundedness_proof(instance :: IPInstance, nonnegative :: Vector{Bool}, i :: Int) :: Vector{Int}

Return a vector `u` in kernel(`instance.A`) proving that the variable of
index `i` is unbounded.
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
        return Int[]
        #error("Unboundedness model should be feasible, status: ",
        #      JuMP.termination_status(model))
    end
    u = Int.(round.(JuMP.value.(vars)))
    return u
end

#
# Permutation-related functions
#

"""
    compute_permutation(bounded :: Vector{Bool}, nonnegative :: Vector{Bool}) :: Tuple{Vector{Int}, Int, Int}

Return a permutation of variables that puts the variables in the order [ bounded ; unbounded and restricted ; unrestricted ], along with the indices of the last bounded variable and the last unbounded but restricted variable.

This operation should be stable with respect to the initial ordering of
variables.

A permutation is represented by a vector perm such that perm[i] = j means
that variable i is sent to j by the permutation.
"""
function compute_permutation(
    bounded::Vector{Bool},
    nonnegative::Vector{Bool}
)::Tuple{Vector{Int},Int,Int}
    #This code is repetitive and inefficient, but I think it is clear
    #For clarity, I'll leave it this way. It is unlikely to become a bottleneck
    @assert length(bounded) == length(nonnegative)
    n = length(bounded)
    permutation = zeros(Int, n)
    #Set bounded variables first after permutation
    n_bounded = 0
    for i in 1:n
        if bounded[i] && nonnegative[i]
            n_bounded += 1
            permutation[n_bounded] = i
        end
    end
    #Next, we need to set non-negative, unbounded variables
    n_nonnegative = 0
    if n_bounded < n #If every variable is bounded, skip this
        for i in 1:n
            if !bounded[i] && nonnegative[i]
                n_nonnegative += 1
                permutation[n_bounded + n_nonnegative] = i
            end
        end
    end
    #Set unrestricted variables last in the permutation
    n_unrestricted = 0
    if n_bounded + n_nonnegative < n
        for i in 1:n
            if !nonnegative[i]
                n_unrestricted += 1
                permutation[n_bounded + n_nonnegative + n_unrestricted] = i
            end
        end
    end
    @assert n_bounded + n_nonnegative + n_unrestricted == n
    bounded_end = n_bounded
    nonnegative_end = n_bounded + n_nonnegative
    return permutation, bounded_end, nonnegative_end
end

"""
    apply_permutation(vector_set :: Vector{Vector{Int}}, permutation :: Vector{Int}) :: Vector{Vector{Int}}

Apply `permutation` to each vector in `vector_set`.
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
    original_variable_order(vector_set :: Vector{Vector{Int}}, instance :: IPInstance) :: Vector{Vector{Int}}

Invert the variable permutation of `instance` over `vector_set`, returning to
the original problem's variable order.

This is useful to give users output in the same variable order they input.
"""
function original_variable_order(
    vector_set :: Vector{Vector{Int}},
    instance :: IPInstance
) :: Vector{Vector{Int}}
    return apply_permutation(vector_set, instance.inverse_permutation)
end

#
# Generating some IPInstances
#

"""
    random_instance(m :: Int, n :: Int) :: IPInstance

Return a random feasible IPInstance with `m` constraints and `n` variables.
"""
function random_ipinstance(
    m :: Int,
    n :: Int
) :: IPInstance
    instance = nothing
    feasible = false
    bounded = false
    while !feasible || !bounded
        #Build random instance in these parameters
        A = rand(-5:5, m, n)
        b = rand(5:20, m)
        C = rand(-10:-1, 1, n)
        u = Union{Int, Nothing}[]
        for _ in 1:n
            push!(u, nothing)
        end
        instance = IPInstance(A, b, C, u, invert_objective=false)
        #Check feasibility
        model, _, _ = SolverTools.feasibility_model(
            instance.A, instance.b, instance.u, nonnegative_vars(instance), Int
        )
        feasible = SolverTools.is_feasible(model)
        SolverTools.set_jump_objective!(model, :Min, vec(instance.C))
        bounded = SolverTools.is_bounded(model)
    end
    return instance
end

end
