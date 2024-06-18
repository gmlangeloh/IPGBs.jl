module IPInstances

export IPInstance, nonnegative_vars, is_bounded, unboundedness_proof, update_objective!, nonnegativity_relaxation, group_relaxation, truncation_weight, projection, project_vector, unbounded_variables, is_feasible_solution, add_constraint, lattice_basis_projection, in_kernel,
apply_permutation, fiber_solution, lattice_basis
export has_slacks, nonnegative_data_only, multiobjective_from_file, num_objectives

using LinearAlgebra

#TODO: Finish removing dependency on AbstractAlgebra. Use only MatrixTools
using AbstractAlgebra
using GLPK
using JuMP

using IPGBs
using IPGBs.GBTools
using IPGBs.MatrixTools
using IPGBs.SolverTools

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
    #Guarantee that A has full row rank
    m, n = size(A)
    if rank(A) < m
        A, b = li_rows(A, b)
    end
    if !apply_normalization
        return A, b, C, u, nonnegative
    end
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
Represents an instance of a problem

min C * x

s.t. A * x = b

0 <= x <= u

x in ZZ^n

The instance is stored in normalized form, with permuted variables so that
the variables appear in the following order: bounded, non-negative but unbounded, unrestricted.
"""
mutable struct IPInstance
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
    binaries :: Vector{Bool}

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
    optimizer :: DataType

    #Lattice-related information
    lattice_basis :: Generic.MatSpaceElem{Int} #Row basis
    rank :: Int
    fiber_solution :: Vector{Int} #v such that Av = b. Not necessarily non-negative.
    computed_lattice_info :: Bool
    originally_bounded :: Vector{Bool}

    #TODO: put a parameter to determine whether it is minimization or not
    function IPInstance(
        A::Array{Int,2},
        b::Vector{Int},
        C::Array{T,2},
        u::Vector{<:Union{Int,Nothing}},
        nonnegative::Union{Nothing,Vector{Bool}} = nothing;
        apply_normalization::Bool = true,
        invert_objective::Bool = true,
        optimizer :: DataType = GLPK.Optimizer
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
        model, model_vars, model_cons = SolverTools.relaxation_model(
            A, b, C, u, nonnegative, optimizer=optimizer
        )
        #Compute a permutation of variables of the given instance such that
        #vars appear in order: bounded, non-negative, unrestricted
        bounded = SolverTools.bounded_variables(A, nonnegative, optimizer=optimizer)
        permutation, bounded_end, nonnegative_end = compute_permutation(bounded, nonnegative)
        inverse_perm = invperm(permutation)
        binaries = [ u[i] == 1 for i in 1:length(u)]
        #Permute columns of problem data
        A = A[:, permutation]
        C = C[:, permutation]
        u = u[permutation]
        nonnegative = nonnegative[permutation]
        #Update the JuMP model with the permutation info
        model, model_vars, model_cons = SolverTools.relaxation_model(
            A, b, C, u, nonnegative, optimizer=optimizer
        )
        #Checks feasibility of the linear relaxation
        @assert SolverTools.is_feasible(model)
        #Compute boundedness of variables using the model
        SolverTools.set_jump_objective!(model, :Min, C[1, :])
        #Initialize lattice data as empty
        basis = AbstractAlgebra.matrix(MatrixTools.AlgebraInt, zeros(Int, m, n))
        rnk = 0
        fiber_sol = Int[]
        #Create the normalized instance
        new(A, b, C, u,
            bounded_end, nonnegative_end, permutation, inverse_perm,
            binaries, m, n, new_m, new_n, true,
            model, model_vars, model_cons, optimizer,
            basis, rnk, fiber_sol, false, bounded
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

function in_kernel(v :: Vector{Int}, instance :: IPInstance)
    return iszero(instance.A * v)
end

function in_kernel(vs :: Vector{Vector{Int}}, instance :: IPInstance)
    return all(v -> in_kernel(v, instance), vs)
end

function multiobjective_from_file(filename :: String)
    parse_line(l) = parse.(Int, split(l))
    f = open(filename, "r")
    line_index = 1
    m, n, p = (0, 0, 0)
    objective_rows = Vector{Int}[]
    constraint_rows = Vector{Int}[]
    rhs = Int[]
    for line in eachline(f)
        data = parse_line(line)
        if line_index == 1
            m = data[1]
            n = data[2]
            p = data[3]
        elseif line_index <= p + 1
            push!(objective_rows, copy(data))
        elseif line_index <= p + m + 1
            push!(constraint_rows, copy(data))
        elseif line_index <= p + m + 2
            rhs = copy(data)
        end
        line_index += 1
    end
    close(f)
    A = foldl(vcat, map(row -> reshape(row, 1, length(row)), constraint_rows))
    C = foldl(vcat, map(row -> reshape(row, 1, length(row)), objective_rows))
    b = rhs
    return IPInstance(A, b, C, fill(nothing, n), apply_normalization=false)
end

function num_objectives(instance :: IPInstance)
    return size(instance.C, 1)
end

function new_first_objective(
    instance :: IPInstance,
    first_objective :: Int
)
    obj_permutation = collect(1:num_objectives(instance))
    obj_permutation[first_objective] = 1
    obj_permutation[1] = first_objective
    return instance.C[obj_permutation, :]
end

"""
    initial_jump_model(
    instance :: IPInstance,
    first_objective :: Int = 1
)

Return a JuMP Model corresponding to the IPInstance.

The JuMP Model keeps all objectives in its objective function matrix.
"""
function jump_model(
    instance :: IPInstance,
    first_objective :: Int = 1;
    optimizer :: DataType = GLPK.Optimizer
)
    C = new_first_objective(instance, first_objective)
    return SolverTools.jump_model(
        instance.A, instance.b, C, instance.u,
        nonnegative_vars(instance), Int, optimizer=optimizer
    )
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

function compute_lattice_info!(instance :: IPInstance)
    if instance.computed_lattice_info
        return
    end
    basis, rnk = hnf_lattice_basis(instance.A)
    fiber_sol = MatrixTools.solve(instance.A, instance.b)
    instance.lattice_basis = basis
    instance.rank = rnk
    instance.fiber_solution = fiber_sol
    instance.computed_lattice_info = true
end

function fiber_solution(instance :: IPInstance)
    compute_lattice_info!(instance)
    return instance.fiber_solution
end

function lattice_basis(instance :: IPInstance)
    compute_lattice_info!(instance)
    return instance.lattice_basis
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

function extract_objective(
    objective_list :: Vector,
    x :: Vector{JuMP.VariableRef}
)
    #Make a matrix whose rows are the coefficients of each objective function
    rows = Matrix{Int}[]
    for obj in objective_list
        c = extract_objective(obj, x)
        push!(rows, c)
    end
    C = foldl(vcat, rows)
    return C
end

function extract_objective(obj :: JuMP.AffExpr, x :: Vector{JuMP.VariableRef})
    c = zeros(Int, 1, length(x))
    for j in 1:length(x)
        if !haskey(obj.terms, x[j])
            continue
        end
        coef = obj.terms[x[j]]
        c[1, j] = round(Int, coef)
    end
    return c
end

function extract_objective(obj::JuMP.VariableRef, x::Vector{JuMP.VariableRef})
    # JuMP represents the objective function as a single VariableRef when
    #possible.
    # This case needs to be treated separately here.
    c = zeros(Int, 1, length(x))
    index = obj.index.value
    c[1, index] = 1
    return c
end

function extract_objective(
    model :: JuMP.Model,
    x :: Vector{JuMP.VariableRef}
) :: Matrix{Int}
    c = extract_objective(objective_function(model), x)
    if objective_sense(model) == MOI.MAX_SENSE
        c = -c
    end
    #If the objective sense is minimization, no normalization of the
    #objective function coefficients is needed.
    return c
end

function IPInstance(
    model::JuMP.Model;
    infer_binary :: Bool = true,
    optimizer :: DataType = GLPK.Optimizer
)
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
                    if !infer_binary
                        continue
                    end
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
    original_rows = length(rows)
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
    #Make sure the original constraints are full row rank
    #And thus that the whole matrix is full row rank. The main reason for doing this
    #with the original constraints is numerical stability.
    frrA, frrB = li_rows(A[1:original_rows, :], b[1:original_rows])
    A = [frrA; A[(original_rows+1):end, :]]
    b = [frrB; b[(original_rows+1):end]]
    #Build upper bound vector
    u :: Vector{Union{Int, Nothing}} = fill(nothing, size(A, 2))
    for (var, ub) in upper_bounds
        u[var] = ub
    end
    #Extend c to the slack variables
    num_slacks = size(A, 2) - size(c, 2)
    Zs = zeros(Int, size(c, 1), num_slacks)
    c = hcat(c, Zs)
    #Build the IPInstance object. The matrices here are already normalized,
    #so no additional normalization is necessary
    return IPInstance(A, b, c, u, apply_normalization=false, optimizer=optimizer)
end

function IPInstance(
    path :: String;
    infer_binary :: Bool = true,
    optimizer :: DataType = GLPK.Optimizer
)
    model = read_from_file(path)
    return IPInstance(model, infer_binary=infer_binary, optimizer=optimizer)
end

"""
    solve(instance :: IPInstance)

    Return the optimal solution to this IPInstance. This solution is computed
    using a traditional IP solver.
"""
function solve(instance :: IPInstance) :: Tuple{Vector{Int}, Int, TerminationStatusCode}
    return SolverTools.solve(
        instance.A, instance.b, instance.C, instance.u,
        nonnegative_variables(instance), Int, optimizer=instance.optimizer
    )
end

function is_feasible_solution(
    instance :: IPInstance,
    solution :: Vector{Int},
    permutation :: Vector{Int} = collect(1:instance.n)
) :: Bool
    @assert length(solution) == instance.n
    perm_solution = solution[permutation]
    return instance.A * perm_solution == instance.b &&
        all(perm_solution[1:instance.nonnegative_end] .>= 0)
end

function guess_initial_solution(
    instance :: IPInstance
) :: Vector{Int}
    A = instance.A
    b = instance.b
    m, n = size(A)
    solution = zeros(Int, n)
    #If the last m columns of A form an identity matrix, they are slack variables
    #This means we can set all other variables to 0 and the slack variables to b,
    #at least if b is non-negative
    if GBTools.ends_with_slacks(A)
        if any(b[i] < 0 for i in 1:m)
            throw(ArgumentError("Cannot guess initial solution for this instance"))
        end
        solution[(n-m+1):n] = instance.b
        return solution
    end
    if n % 2 == 0
        #Guess that half of the variables are binary slacks. Then if we have a guess for
        #the other half, we can deduce the binary slacks and get a full solution
        half = n ÷ 2
        if isinteger(sqrt(half)) #Guess: assignment problem
            #Permutation solution.
            k = Int(sqrt(half))
            for i in 1:k
                solution[i + (i - 1) * k] = 1
            end
            #Complete this partial solution with binary slacks
            for i in half+1:n
                solution[i] = 1 - solution[i - half]
            end
            return solution
        end
    end
    throw(ArgumentError("Cannot guess initial solution for this instance"))
end

function extend_feasible_solution(instance :: IPInstance, solution :: Vector{Int})
    return MatrixTools.lift_partial_solution(solution, instance.b, instance.A)
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

has_slacks(instance :: IPInstance) = GBTools.has_slacks(instance.A)

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

function unbounded_variables(instance :: IPInstance) :: Vector{Bool}
    return [ i > instance.bounded_end && i <= instance.nonnegative_end for i in 1:instance.n]
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
    # Update the objective function as well.
    # For the definition, see Thomas - The Structure of Group Relaxations
    #If we are relaxing a full rank set of non-negative variables, we can
    #update the objective function. Otherwise, we are in a strictly extended
    #group relaxation, and we can keep the objective function of the previous
    #relaxation.
    new_C = instance.C
    if length(nonnegative) - count(nonnegative) == instance.m
        relaxed = [!nonnegative[i] for i in 1:instance.n]
        inv_As = inv(instance.A[:, relaxed])
        Ans = instance.A[:, nonnegative]
        cs = instance.C[1, relaxed]
        cns = instance.C[1, nonnegative]
        new_obj = cns' - cs' * inv_As * Ans
        new_C = zeros(Float64, size(instance.C))
        j = 1
        for i in 1:length(nonnegative)
            if nonnegative[i]
                new_C[1, i] = new_obj[j]
                j += 1
            end
        end
    end
    return IPInstance(
        instance.A, instance.b, new_C, instance.u, nonnegative,
        apply_normalization=false,
        invert_objective=false
    )
end

function complement(
    s :: Vector{Int},
    n :: Int
) :: Vector{Int}
    comp = Int[]
    for i in 1:n
        if i in s
            continue
        end
        push!(comp, i)
    end
    return comp
end

function projection(
    instance :: IPInstance,
    away_from :: Vector{Int}
) :: IPInstance
    #Only project away from variables with relaxed non-negativity
    @assert all(s > instance.nonnegative_end for s in away_from)
    #Find variables to project onto
    onto = complement(away_from, instance.n)
    #Project onto the variables in onto
    #TODO: Fix this. Just projecting the constraint matrix doesn't work.
    return IPInstance(
        instance.A[:, onto], instance.b, instance.C[:, onto], instance.u[onto],
        apply_normalization=false,
        invert_objective=false
    )
end

function project_vector(
    v :: Vector{Int},
    away_from :: Vector{Int}
)
    onto = complement(away_from, length(v))
    return v[onto]
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
    var_basis = SolverTools.optimal_basis!(instance.model, instance.model_vars)
    nonbasics = [ !variable for variable in var_basis ]
    @assert count(var_basis) == instance.m
    #Keep only the nonnegativity constraints on the non-basic variables
    return nonnegativity_relaxation(instance, nonbasics)
end

function lattice_basis_projection(
    instance :: IPInstance,
    var_selection :: Symbol = :Any
)
    if var_selection == :Any
        li_cols = Int[]
        sigma = Int[] #Complement of li_cols
        #Greedy approach: pick a column and then check for linear independence.
        #For simplicity, I can check this by looking at whether the rank increased
        L = lattice_basis(instance)
        j = 1
        while j <= instance.n
            push!(li_cols, j)
            #Check linear independence
            li_basis = L[:, li_cols]
            if rank(li_basis) < length(li_cols)
                #If not linearly independent, remove the last column
                pop!(li_cols)
                push!(sigma, j)
            end
            j += 1
        end
    elseif var_selection == :SimplexBasis
        #Create some LI set for the lattice from the non-basic
        #variables of the optimal solution to the linear relaxation
        var_basis = SolverTools.optimal_basis!(instance.model, instance.model_vars)
        sigma = [ i for i in eachindex(var_basis) if var_basis[i] ]
        new_rank = rank(instance.A[:, sigma])
        old_rank = new_rank
        #Something weird happened numerically with the Simplex basis. We may need
        #to delete a part of it.
        while new_rank < count(var_basis)
            #Find some element of var_basis that can be deleted while keeping the same rank
            for j in 1:length(var_basis)
                if !var_basis[j]
                    continue
                end
                var_basis[j] = false
                new_sigma = [ i for i in eachindex(var_basis) if var_basis[i] ]
                old_rank = new_rank
                new_rank = rank(instance.A[:, new_sigma])
                if new_rank != old_rank
                    var_basis[j] = true
                else #We found a variable to delete
                    break
                end
            end
        end
        #At worst, we may be missing one variable in var_basis to get to the full rank of A
        #That happens when the objective function is linearly dependent with some constraint
        #In that case, we can just add some new LI variable not in var_basis
        old_rank = new_rank
        j = 1
        while j <= instance.n && old_rank < instance.m
            push!(sigma, j)
            new_rank = rank(instance.A[:, sigma])
            if !var_basis[j] && new_rank > old_rank
                var_basis[j] = true
                old_rank = new_rank
            else
                pop!(sigma)
            end
            j += 1
        end
        sort!(sigma)
        li_cols = [ i for i in eachindex(var_basis) if !var_basis[i] ]
    else
        error("Unknown variable selection method: $var_selection")
    end
    basis = lattice_basis(instance)[:, li_cols]
    return basis_to_uhnf(basis), basis, sigma
end

function linear_relaxation(instance :: IPInstance)
    optimize!(instance.model)
    return objective_value(instance.model)
end

function truncation_weight(
    instance :: IPInstance
) :: Tuple{Vector{Float64}, Float64}
    A = Array(Int.(lattice_basis(instance)))
    b = fiber_solution(instance)
    unbounded = map(x -> !x, instance.originally_bounded)
    return SolverTools.optimal_weight_vector(A, b, unbounded, optimizer=instance.optimizer)
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

function update_objective!(
    instance :: IPInstance,
    j :: Int,
    sigma :: Vector{Int}
)
    c = SolverTools.bounded_objective(instance.A, j, sigma, optimizer=instance.optimizer)
    #We take the negative here to normalize the problem to minimization form
    instance.C[1, :] = c
    #To ease debugging: check whether c has the specified properties, that is,
    #1. c[sigma] == 0 and
    #2. c' * u == -u[i] for all u in the lattice basis
    @assert all(IPGBs.is_approx_zero.(c[sigma]))
    for i in 1:size(lattice_basis(instance), 1)
        u = reshape(Array(lattice_basis(instance)[i, :]), instance.n)
        @assert abs(c' * u + u[j]) < 1e-6
    end
    #Update the LP model as well
    SolverTools.set_jump_objective!(instance.model, :Min, instance.C[1, :])
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
    unboundedness_proof(
    instance :: IPInstance,
    i :: Int
) :: Vector{Int}

Return a vector `u` in kernel(`instance.A`) proving that the variable of
index `i` is unbounded.
"""
function unboundedness_proof(
    instance :: IPInstance,
    i :: Int
) :: Vector{Int}
    model, vars, _ = SolverTools.unboundedness_ip_model(
        instance.A, nonnegative_variables(instance), i, optimizer=instance.optimizer
    )
    JuMP.optimize!(model)
    if JuMP.termination_status(model) != JuMP.MOI.OPTIMAL
        return Int[]
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

function add_constraint(
    instance :: IPInstance,
    constraint :: Vector{Int},
    rhs :: Int,
    sense :: Symbol = :EQ
) :: IPInstance
    #Take the constraint matrix, RHS and costs, then add a new slack variable
    #for the new constraint if it's an inequality
    A = instance.A
    b = instance.b
    C = instance.C
    u = instance.u
    new_A = [ A ; constraint' ]
    new_b = [ b ; rhs ]
    new_C = C
    new_u = u
    if sense == :EQ
        #No need to add a slack variable
    elseif sense == :LT || sense == :GT
        new_A = [ new_A zeros(Int, size(new_A, 1), 1) ]
        new_A[end, end] = sense == :LT ? 1 : -1
        new_C = [ C zeros(Float64, size(C, 1), 1) ]
        new_u = [ u ; nothing ] #TODO Check whether I want an explicit upper bound
    else
        error("Unknown constraint sense: $sense")
    end
    return IPInstance(new_A, new_b, new_C, new_u, apply_normalization=false)
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
    #TODO: Put this somewhere else!!!
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
            instance.A, instance.b, instance.u, nonnegative_vars(instance), Int,
            optimizer=instance.optimizer
        )
        feasible = SolverTools.is_feasible(model)
        SolverTools.set_jump_objective!(model, :Min, vec(instance.C))
        bounded = SolverTools.is_bounded(model)
    end
    return instance
end

end
