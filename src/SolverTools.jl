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
Generates a feasibility checking model for

Ax = b
0 <= x <= u

where x is either an integer variable vector, if `var_type` == Int,
or a real variable vector otherwise.

Returns the JuMP model along with vectors with references to its variables
and constraints.

TODO should I do something special for the binary case?
"""
function feasibility_model(
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int},
    var_type :: DataType
) :: Tuple{JuMP.Model, Vector{VariableRef}, Vector{ConstraintRef}}
    m, n = size(A)
    model = Model(GENERAL_SOLVER.Optimizer)
    set_silent(model)
    if var_type == Int #use the original IP, not the linear relaxation
        @variable(model, x[1:n] >= 0, Int)
    else #var_type is Real / linear relaxation is used
        @variable(model, x[1:n] >= 0)
    end
    set_upper_bound.(x, u)
    #Set constraints, keeping references in a vector for later changes
    constraints = []
    for i in 1:m
        ai = A[i, :]
        con = @constraint(model, ai' * x == b[i])
        push!(constraints, con)
    end
    @objective(model, Max, 0)
    return model, x, constraints
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
