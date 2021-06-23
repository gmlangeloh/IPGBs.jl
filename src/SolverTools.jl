"""
Includes various LP and IP functions using external solvers and JuMP.
"""
module SolverTools

using JuMP
using Clp
using CPLEX #TODO have to check license stuff for distribution later...
#What do I do to select another solver in case CPLEX is not installed?

const LP_SOLVER = Clp
const GENERAL_SOLVER = CPLEX

"""
Computes a strictly positive vector in the row span of `A` using
linear programming. Assumes Ax = b is feasible, and
max x
s.t. Ax = b
x >= 0
is bounded. Given these conditions, the dual variables of the
constraints of the above LP give a positive row span vector.
"""
function positive_row_span(
    A :: Array{Int, 2},
    b :: Vector{Int}
) :: Vector{Float64}
    m, n = size(A)
    model = Model(GENERAL_SOLVER.Optimizer)
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

"""
Returns a JuMP model (alongside references to its variables and constraints)
for the problem

min c[1, :] * x
s.t. Ax = b
x_i >= 0 for all i s.t. nonnegative[i] == true

TODO Should I do anything special for the binary case?
"""
function jump_model(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Float64},
    u :: Vector{Int},
    nonnegative :: Vector{Bool},
    var_type :: DataType
) :: Tuple{JuMP.Model, Vector{VariableRef}, Vector{ConstraintRef}}
    m, n = size(A)
    model = Model(GENERAL_SOLVER.Optimizer)
    set_silent(model)
    if var_type == Int #use the original IP, not the linear relaxation
        @variable(model, x[1:n], Int)
    else #var_type is Real / linear relaxation is used
        @variable(model, x[1:n])
    end
    set_upper_bound.(x, u)
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

"""
Returns a linear relaxation model of the problem
min C[1, :]' * x
s.t. Ax = b
x >= 0
"""
function relaxation_model(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Float64},
    u :: Vector{Int},
    nonnegative :: Vector{Bool}
) :: Tuple{JuMP.Model, Vector{VariableRef}, Vector{ConstraintRef}}
    return jump_model(A, b, C, u, nonnegative, Real)
end

"""
Searches for integer u such that
Au = 0
u_{not sigma} >= 0
u_i > 0

TODO this could also be done with LP as follows:
"Assume all data is rational. Then, the polyhedron is rational, so the optimum
must be rational. Multiply by a large enough integer..."
Implement it this way later!
"""
function unboundedness_ip_model(
    A :: Array{Int, 2},
    nonnegative :: Vector{Bool},
    i :: Int
) :: Tuple{JuMP.Model, Vector{VariableRef}, Vector{ConstraintRef}}
    #Get model with 0 in RHS and objective function
    m, n = size(A)
    b = zeros(Int, m)
    C = zeros(Int, 1, n)
    u = [typemax(Int) for _ in 1:n]
    model, vars, constrs = jump_model(A, b, C, u, nonnegative, Int)
    @constraint(model, vars[i] >= 1)
    return model, vars, constrs
end

"""
Returns true iff maximizing x_i in model is bounded.

It is assumed that the given model is (primal) feasible.
"""
function is_bounded(
    i :: Int,
    model :: JuMP.Model,
    x :: Vector{VariableRef}
) :: Bool
    @objective(model, Max, x[i])
    optimize!(model)
    #Note that this only guarantees unboundedness if the model
    #is feasible.
    if termination_status(model) != MOI.DUAL_INFEASIBLE
        return true
    end
    return false
end

"""
Generates a feasibility checking model for

Ax = b
0 <= x <= u

where x is either an integer variable vector, if `var_type` == Int,
or a real variable vector otherwise.

Returns the JuMP model along with vectors with references to its variables
and constraints.
"""
function feasibility_model(
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int},
    nonnegative :: Vector{Bool},
    var_type :: DataType
) :: Tuple{JuMP.Model, Vector{VariableRef}, Vector{ConstraintRef}}
    feasibility_obj = zeros(Float64, 1, size(A, 2))
    return jump_model(A, b, feasibility_obj, u, nonnegative, var_type)
end

"""
Changes the RHS of `constraints` to b - A * v.
"""
function update_feasibility_model_rhs(
    constraints :: Vector{ConstraintRef},
    A :: Array{Int, 2},
    b :: Vector{Int},
    v :: T
) where {T <: AbstractVector{Int}}
    delta = A * v
    #new_rhs = (b - delta)[1:length(constraints)]
    set_normalized_rhs.(constraints, b - delta)
end

"""
Returns true iff the IP/LP given by `model` is feasible.
"""
function is_feasible(
    model :: JuMP.Model
) :: Bool
    optimize!(model)
    if termination_status(model) == MOI.INFEASIBLE
        return false
    end
    return true
end

end
