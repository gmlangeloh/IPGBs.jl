"""
Includes various LP and IP functions using external solvers and JuMP.
"""
module SolverTools

using JuMP
using GLPK

using IPGBs

const SOLVER = GLPK.Optimizer

"""
    cone_element(
    rays :: Vector{Vector{T}}
) :: Vector{Float64} where {T <: Real}

Return a point in the cone spanned by the rays in `rays` if it exists.

The point is computed via linear programming.
"""
function cone_element(
    rays :: Vector{Vector{T}};
    optimizer = SOLVER
) :: Vector{Float64} where {T <: Real}
    if isempty(rays)
        return Float64[]
    end
    n = length(rays[1])
    model = Model(optimizer)
    @variable(model, x[1:n] >= 0)
    for ray in rays
        @constraint(model, ray' * x >= 0)
    end
    @constraint(model, sum(x) == 1)
    @objective(model, Min, sum(x))
    set_silent(model)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return value.(x)
    end
    #The cone is empty
    return Float64[]
end

"""
    optimal_basis!(model :: JuMP.Model) :: Vector{Bool}

Return a boolean vector with value true at index i iff the i-th variable is
basic at the optimal solution of `model`.

This function calls `optimize!` on `model` in order to be self-contained.
"""
function optimal_basis!(
    model::JuMP.Model,
    x :: Vector{JuMP.VariableRef}
)::Vector{Bool}
    optimize!(model)
    #Reoptimize with the objective value set as a constraint.
    #The new objective function will be a sum of the values of the original variables
    #This guarantees an optimal basis involving only these variables, and not slacks
    #corresponding to the constraints.
    val = objective_value(model)
    old_obj = objective_function(model)
    obj_constr = @constraint(model, old_obj == val)
    @objective(model, Max, sum(x))
    optimize!(model)
    #Find basic original variables for the problem
    n = length(x)
    var_basis = fill(false, n)
    found_zero = MOI.get(model, MOI.ConstraintBasisStatus(), obj_constr) == MOI.BASIC
    for j in 1:n
        status = MOI.get(model, MOI.VariableBasisStatus(), x[j])
        #It is necessary to eliminate a basic variable with value 0, because we added
        #a new constraint to the model, so the optimal basis will contain an extra
        #basic variable. Such a variable should always exist.
        if status == MOI.BASIC && isapprox(value(x[j]), 0.0, atol=1e-8) && !found_zero
            found_zero = true
            continue
        end
        if status == MOI.BASIC
            var_basis[j] = true
        end
    end
    #Go back to the original model
    @objective(model, Min, old_obj)
    delete(model, obj_constr)
    return var_basis
end

"""
    is_degenerate(model :: JuMP.Model, vars :: Vector{VariableRef}, constraints :: Vector{ConstraintRef})

Return true iff the optimal basis of `model` is degenerate, i.e. if any of its
basic variables has value 0.
"""
function is_degenerate(
    model::JuMP.Model,
    vars::Vector{JuMP.VariableRef},
    constraints::Vector{JuMP.ConstraintRef}
)::Bool
    for i in 1:length(vars)
        if MOI.get(model, MOI.VariableBasisStatus(), vars[i]) == MOI.BASIC
            if IPGBs.is_approx_zero(value(vars[i]))
                return true
            end
        end
    end
    #TODO: Check at some point whether this really corresponds to the slacks
    #of the constraints!
    for j in 1:length(constraints)
        if MOI.get(model, MOI.ConstraintBasisStatus(), constraints[j]) == MOI.BASIC
            if IPGBs.is_approx_zero(value(constraints[j]))
                return true
            end
        end
    end
    return false
end

"""
    optimal_row_span(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{T},
    sense::Symbol = :Min
)::Union{Vector{Float64},Nothing} where {T<:Real}

Compute a vector in the row span of `A` given by `yA` where `y` is the optimal
solution to the dual of the LP given by `A, b, c`.

Assumes Ax = b is feasible, and max x s.t. Ax = b, x >= 0 is bounded.
Given these conditions, the dual variables of the constraints of the above
LP give a positive row span vector.
"""
function optimal_row_span(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{T},
    sense::Symbol = :Min;
    optimizer = SOLVER
) where {T<:Real}
    m, n = size(A)
    @assert(n == size(C, 2))
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n] >= 0)
    constraints = []
    for i in 1:m
        ai = A[i, :]
        con = @constraint(model, ai' * x == b[i])
        push!(constraints, con)
    end
    if sense == :Max
        @objective(model, Max, (C*x)[1])
    else
        @objective(model, Min, (C*x)[1])
    end
    optimize!(model)
    if !has_duals(model)
        return nothing
    end
    return A' * shadow_price.(constraints)
end

"""
    positive_row_span(A :: Matrix{Int}, b :: Vector{Int})

Compute a strictly positive vector in the row span of `A` using
linear programming.
"""
function positive_row_span(
    A::Matrix{Int},
    b::Vector{Int};
    optimizer=SOLVER
)::Union{Vector{Float64},Nothing}
    obj = ones(Float64, 1, size(A, 2))
    return optimal_row_span(A, b, obj, :Max, optimizer=optimizer)
end

"""
    jump_model(A :: Matrix{Int}, b :: Vector{Int}, C :: Array{Float64}, u :: Vector{Union{Int, Nothing}}, nonnegative :: Vector{Bool}, var_type :: DataType)

Return a JuMP model (alongside references to its variables and constraints)
for the given IP or LP problem.

The optimization problem considered is
min c[1, :] * x s.t. Ax = b, x_i >= 0 for all i s.t. nonnegative[i] == true

The variables are integer if `var_type == Int` or real otherwise.

TODO: Should I do anything special for the binary case?
"""
function jump_model(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{Float64},
    u::Vector{<: Union{Int, Nothing}},
    nonnegative::Vector{Bool},
    var_type::DataType;
    optimizer = SOLVER
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    m, n = size(A)
    model = Model(optimizer)
    set_silent(model)
    if var_type == Int #use the original IP, not the linear relaxation
        @variable(model, x[1:n], Int)
    else #var_type is Real / linear relaxation is used
        @variable(model, x[1:n])
    end
    #TODO: Should I set these upper bounds? My model already includes them
    #when necessary as additional constraints!
    #for i in 1:n
    #    if !isnothing(u[i])
    #        set_upper_bound(x[i], u[i])
    #    end
    #end
    #Set non-negativity constraints for the relevant variables
    for i in 1:n
        if nonnegative[i]
            set_lower_bound(x[i], 0)
        end
    end
    #Set constraints, keeping references in a vector for later changes
    constraints = []
    for i in 1:m
        ai = A[i, :]
        con = @constraint(model, ai' * x == b[i])
        push!(constraints, con)
    end
    @objective(model, Min, C[1, :]' * x)
    return model, x, constraints
end

function solve(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{Float64},
    u::Vector{<: Union{Int, Nothing}},
    nonnegative::Vector{Bool},
    var_type::DataType;
    optimizer=SOLVER
):: Tuple{Vector{Int}, Int, TerminationStatusCode}
    model, x, _ = jump_model(A, b, C, u, nonnegative, var_type, optimizer=optimizer)
    optimize!(model)
    if !has_values(model)
        return zeros(Int, length(x)), 0, termination_status(model)
    end
    solution = round.(Int, value.(x))
    val = round(Int, objective_value(model))
    return solution, val, termination_status(model)
end

"""
    relaxation_model(A :: Matrix{Int}, b :: Vector{Int}, C :: Array{Float64}, u :: Vector{Union{Int, Nothing}}, nonnegative :: Vector{Bool})

Return a linear relaxation model of the given IP.

The optimization problem considered is of the form
min C[1, :]' * x s.t. Ax = b, x >= 0.
"""
function relaxation_model(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{Float64},
    u::Vector{<: Union{Int, Nothing}},
    nonnegative::Vector{Bool};
    optimizer = SOLVER
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    return jump_model(A, b, C, u, nonnegative, Real, optimizer=optimizer)
end

"""
    unboundedness_ip_model(A :: Matrix{Int}, nonnegative :: Vector{Bool}, i :: Int)

Return an IP that checks for the existence of an integer vector `u` such that
`Au = 0`, `u[nonnegative] >= 0`, `u[i] > 0`.

TODO: this could also be done with LP as follows:
"Assume all data is rational. Then, the polyhedron is rational, so the optimum
must be rational. Multiply by a large enough integer..."
Implement it this way later!
"""
function unboundedness_ip_model(
    A::Array{Int,2},
    nonnegative::Vector{Bool},
    i::Int;
    optimizer = SOLVER
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    #Get model with 0 in RHS and objective function
    m, n = size(A)
    b = zeros(Int, m)
    C = zeros(Float64, 1, n)
    u = [nothing for _ in 1:n]
    model, vars, constrs = jump_model(A, b, C, u, nonnegative, Int, optimizer=optimizer)
    @constraint(model, vars[i] >= 1)
    return model, vars, constrs
end

function bounded_variables(
    A :: Matrix{Int},
    nonnegative :: Vector{Bool};
    optimizer = SOLVER
)
    m, n = size(A)
    bounded = Bool[]
    #Find some feasible RHS. This does not depend on the variable
    #we are trying to bound.
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n])
    @variable(model, y[1:m])
    for k in 1:n
        if nonnegative[k]
            @constraint(model, x[k] >= 0)
        end
    end
    @objective(model, Min, 0)
    for k in 1:m
        @constraint(model, A[k, :]' * x == y[k])
    end
    optimize!(model)
    b = value.(y)
    #Check boundedness for the RHS found above
    #Bounded for some feasible RHS = bounded for all feasible RHS
    #This has to be done for each variable separately, but we can reuse
    #the JuMP LP model for efficiency.
    opt_model = Model(SOLVER)
    set_silent(opt_model)
    @variable(opt_model, z[1:n])
    for k in 1:n
        if nonnegative[k]
            @constraint(opt_model, z[k] >= 0)
        end
    end
    for k in 1:m
        @constraint(opt_model, A[k, :]' * z == b[k])
    end
    for i in 1:n
        bnd = true
        @objective(opt_model, Max, z[i])
        optimize!(opt_model)
        if termination_status(opt_model) == MOI.DUAL_INFEASIBLE ||
            termination_status(opt_model) == MOI.INFEASIBLE_OR_UNBOUNDED
            bnd = false
        end
        push!(bounded, bnd)
    end
    return bounded
end

"""
    set_jump_objective!(model :: JuMP.Model, direction :: Symbol, c :: Vector{T}, x :: Vector{JuMP.VariableRef}) where {T <: Real}

Updates `model` changing its objective function to `c` and sense to
`direction`, which can be either `:Max` or `:Min`.
"""
function set_jump_objective!(
    model::JuMP.Model,
    direction::Symbol,
    c::Vector{T},
    x :: Vector{JuMP.VariableRef}
) where {T<: Real}
    if direction == :Max
        @objective(model, Max, c' * x)
    else
        @objective(model, Min, c' * x)
    end
end

function set_jump_objective!(
    model::JuMP.Model, direction::Symbol, c::Vector{T}
) where {T<:Real}
    set_jump_objective!(model, direction, c, all_variables(model))
end

"""
    feasibility_model(A :: Matrix{Int}, b :: Vector{Int}, u :: Vector{<: Union{Int, Nothing}}, nonnegative :: Vector{Bool}, var_type :: DataType)

Return a feasibility checking model along with variable and constraint vectors
for Ax = b, 0 <= x <= u, where x is either an integer variable vector, if
`var_type == Int` or a real variable vector otherwise.
"""
function feasibility_model(
    A::Array{Int,2},
    b::Vector{Int},
    u::Vector{<: Union{Int, Nothing}},
    nonnegative::Vector{Bool},
    var_type::DataType;
    optimizer = SOLVER
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    feasibility_obj = zeros(Float64, 1, size(A, 2))
    return jump_model(A, b, feasibility_obj, u, nonnegative, var_type, optimizer=optimizer)
end

"""
    update_feasibility_model_rhs(constraints :: Vector{ConstraintRef}, A :: Matrix{Int}, b :: Vector{Int}, v :: T) where {T <: AbstractVector{Int}}

Change the RHS of `constraints` to `b - A * v`.
"""
function update_feasibility_model_rhs(
    constraints::Vector{ConstraintRef},
    A::Array{Int,2},
    b::Vector{Int},
    v::T
) where {T<:AbstractVector{Int}}
    delta = A * v
    #new_rhs = (b - delta)[1:length(constraints)]
    set_normalized_rhs.(constraints, b - delta)
end

"""
    is_feasible(model :: JuMP.Model)

Return true iff the IP/LP given by `model` is feasible.
"""
function is_feasible(
    model::JuMP.Model
)::Bool
    optimize!(model)
    if termination_status(model) == MOI.INFEASIBLE
        return false
    end
    return true
end

"""
    is_bounded(model :: JuMP.Model)

Return true iff the IP/LP given by `model` is bounded.

`model` is assumed to be feasible.
"""
function is_bounded(
    model :: JuMP.Model
) :: Bool
    optimize!(model)
    if termination_status(model) in [MOI.DUAL_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        return false
    end
    return true
end

function is_bounded_polyhedron(
    A :: Matrix{Int};
    optimizer = SOLVER
) :: Bool
    m, n = size(A)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:n] >= 0)
    @constraint(model, A * x == zeros(Int, m))
    optimize!(model)
    return termination_status(model) == MOI.OPTIMAL
end

"""
    bounded_objective(A :: Matrix{Int}, i :: Int, sigma :: Vector{Int})

Compute the objective function for the bounded case of the project-and-lift
algorithm.

By Farkas Lemma, either
A^σ x^σ + A_σ x_σ = 0 has a solution with b^T x < 0
or
y^T A^σ <= (b^σ)^T
y^T A_σ == b_σ^T has a solution.
In this case of P&L, the latter holds and if y is a solution it follows that
c = b - A^T y satisfies:
c_σ = 0 and c^T u = -u_i for all u in ker(A).
"""
function bounded_objective(A::Matrix{Int}, i::Int, sigma::Vector{Int}; optimizer = SOLVER)
    m, n = size(A)
    model = Model(optimizer)
    set_silent(model)
    @variable(model, x[1:m])
    @objective(model, Max, 0)
    #TODO: this can be done more efficiently, but it's fine like this for now
    Aσ = A[:, sigma]
    not_sigma = [j for j in 1:n if !(j in sigma)]
    ANσ = A[:, not_sigma]
    b = zeros(Int, n)
    b[i] = -1
    for j in eachindex(sigma)
        @constraint(model, Aσ[:, j]' * x == b[sigma[j]])
    end
    for j in eachindex(not_sigma)
        @constraint(model, ANσ[:, j]' * x <= b[not_sigma[j]])
    end
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        c = value.(x)
        return b - transpose(A) * c
    end
    error("This problem was supposed to be feasible.")
end

"""
    optimal_weight_vector(
    A :: Matrix{Int},
    b :: Vector{Int},
    unbounded :: Vector{Bool}
) :: Tuple{Vector{Float64}, Float64}

    Find an extreme ray of the dual cone of Ax = 0 whenever possible.

    See Malkin's thesis, page 83 for a description of the use of this
    idea in Gröbner basis truncation.
"""
function optimal_weight_vector(
    A :: Matrix{Int},
    b :: Vector{Int},
    unbounded :: Vector{Bool};
    optimizer = SOLVER
) :: Tuple{Vector{Float64}, Float64}
    model = Model(optimizer)
    set_silent(model)
    m, n = size(A)
    @assert length(b) == n
    @variable(model, x[1:n] >= 0)
    @objective(model, Min, b' * x)
    #TODO: Understand whether this is necessary or not.
    #for i in 1:n
    #    if unbounded[i]
    #        set_upper_bound(x[i], 0)
    #    end
    #end
    for i in 1:m
        @constraint(model, A[i, :]' * x == 0)
    end
    @constraint(model, sum(x[j] for j in 1:n) == 1)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return value.(x), objective_value(model)
    end
    #If no truncation vector is found, simply return 0 to
    #disable weight truncation.
    return zeros(Float64, n), 0.0
end

end
