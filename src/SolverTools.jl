"""
Includes various LP and IP functions using external solvers and JuMP.
"""
module SolverTools

using JuMP
using Clp
using CPLEX

using IPGBs

const LP_SOLVER = Clp
const GENERAL_SOLVER = CPLEX

"""
    optimal_basis!(model :: JuMP.Model) :: Vector{Bool}

Return a boolean vector with value true at index i iff the i-th variable is
basic at the optimal solution of `model`.

This function calls `optimize!` on `model` in order to be self-contained.
"""
function optimal_basis!(
    model::JuMP.Model
)::Vector{Bool}
    optimize!(model)
    x = all_variables(model)
    n = length(x)
    basis = fill(false, n)
    for i in 1:n
        status = MOI.get(model, MOI.VariableBasisStatus(), x[i])
        if status == MOI.BASIC
            basis[i] = true
        end
    end
    return basis
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
    #TODO Check at some point whether this really corresponds to the slacks
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
    optimal_row_span(A :: Matrix{Int}, b :: Vector{Int}, c :: Array{T}) where {T <: Real}

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
    sense::Symbol = :Min
)::Union{Vector{Float64},Nothing} where {T<:Real}
    m, n = size(A)
    @assert(n == size(C, 2))
    model = Model(GENERAL_SOLVER.Optimizer)
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
    b::Vector{Int}
)::Union{Vector{Float64},Nothing}
    obj = ones(Float64, 1, size(A, 2))
    return optimal_row_span(A, b, obj, :Max)
end

"""
    jump_model(A :: Matrix{Int}, b :: Vector{Int}, C :: Array{Float64}, u :: Vector{Int}, nonnegative :: Vector{Bool}, var_type :: DataType)

Return a JuMP model (alongside references to its variables and constraints)
for the given IP or LP problem.

The optimization problem considered is
min c[1, :] * x s.t. Ax = b, x_i >= 0 for all i s.t. nonnegative[i] == true

The variables are integer if `var_type == Int` or real otherwise.

TODO Should I do anything special for the binary case?
"""
function jump_model(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{Float64},
    u::Vector{Int},
    nonnegative::Vector{Bool},
    var_type::DataType
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    m, n = size(A)
    model = Model(GENERAL_SOLVER.Optimizer)
    set_silent(model)
    if var_type == Int #use the original IP, not the linear relaxation
        @variable(model, x[1:n], Int)
    else #var_type is Real / linear relaxation is used
        @variable(model, x[1:n])
    end
    for i in 1:n
        if u[i] != typemax(Int) #No upper bound necessary when u[i] == max int
            set_upper_bound(x[i], u[i])
        end
    end
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
    relaxation_model(A :: Matrix{Int}, b :: Vector{Int}, C :: Array{Float64}, u :: Vector{Int}, nonnegative :: Vector{Bool})

Return a linear relaxation model of the given IP.

The optimization problem considered is of the form
min C[1, :]' * x s.t. Ax = b, x >= 0.
"""
function relaxation_model(
    A::Matrix{Int},
    b::Vector{Int},
    C::Array{Float64},
    u::Vector{Int},
    nonnegative::Vector{Bool}
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    return jump_model(A, b, C, u, nonnegative, Real)
end

"""
    unboundedness_ip_model(A :: Matrix{Int}, nonnegative :: Vector{Bool}, i :: Int)

Return an IP that checks for the existence of an integer vector `u` such that
`Au = 0`, `u[nonnegative] >= 0`, `u[i] > 0`.

TODO this could also be done with LP as follows:
"Assume all data is rational. Then, the polyhedron is rational, so the optimum
must be rational. Multiply by a large enough integer..."
Implement it this way later!
"""
function unboundedness_ip_model(
    A::Array{Int,2},
    nonnegative::Vector{Bool},
    i::Int
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
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
    is_bounded(i :: Int, model :: JuMP.Model, x :: Vector{JuMP.VariableRef})

Return true iff maximizing x_i in model is bounded (has an optimal solution).

It is assumed that the given model is (primal) feasible.
"""
function is_bounded(
    i::Int,
    model::JuMP.Model,
    x::Vector{VariableRef}
)::Bool
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

set_jump_objective!(model::JuMP.Model, direction::Symbol, c::Vector{T}) where {T<:Real} = set_jump_objective!(model, direction, c, all_variables(model))

"""
    feasibility_model(A :: Matrix{Int}, b :: Vector{Int}, u :: Vector{Int}, nonnegative :: Vector{Bool}, var_type :: DataType)

Return a feasibility checking model along with variable and constraint vectors
for Ax = b, 0 <= x <= u, where x is either an integer variable vector, if
`var_type == Int` or a real variable vector otherwise.
"""
function feasibility_model(
    A::Array{Int,2},
    b::Vector{Int},
    u::Vector{Int},
    nonnegative::Vector{Bool},
    var_type::DataType
)::Tuple{JuMP.Model,Vector{VariableRef},Vector{ConstraintRef}}
    feasibility_obj = zeros(Float64, 1, size(A, 2))
    return jump_model(A, b, feasibility_obj, u, nonnegative, var_type)
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

Return true iff the IP/LP given by `model` is feasible.
"""
function is_bounded(
    model :: JuMP.Model
) :: Bool
    optimize!(model)
    if termination_status(model) == MOI.DUAL_INFEASIBLE
        return false
    end
    return true
end

end
