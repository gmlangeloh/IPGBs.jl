"""
Some useful functions to deal with multi-objective optimization problems.
"""
module MultiObjectiveTools
export check_size_consistency, is_nondominated, is_nondominated_set, contained_by,
is_efficient_set_feasible, remove_dominated!, ideal_point, nadir_bound

using IPGBs
using IPGBs.IPInstances

using JuMP
using GLPK

function check_size_consistency(
    A :: Matrix{Int},
    b :: Vector{Int},
    C :: Matrix{Int},
    x :: Vector{Int}
)
    @assert size(A, 2) == size(C, 2)
    @assert size(A, 2) == length(x)
    @assert size(A, 1) == length(b)
end

"""
Creates a vector of pareto points from a set of pareto points
"""
function pareto_vector(
    pareto :: Set{Vector{Int}}
) :: Vector{Vector{Int}}
    pareto_vec = Vector{Int}[]
    for nondominated in pareto
        push!(pareto_vec, nondominated)
    end
    return pareto_vec
end

"""
Checks whether the solution at index `i` of `pareto` is nondominated with respect
to the rest of `pareto`.

Returns `nothing` in case the solution is nondominated, or an index pointing to
proof that the solution is inefficient, otherwise.
"""
function is_nondominated(
    i :: Int,
    pareto :: Vector{Vector{Int}};
    maximization :: Bool = true
) :: Union{Int, Nothing}
    if isempty(pareto)
        return nothing
    end
    @assert i <= length(pareto)
    n = length(pareto[i])
    for j in 1:length(pareto)
        if i != j
            if maximization
                if all(k -> pareto[i][k] <= pareto[j][k], 1:n)
                    return j
                end
            else
                if all(k -> pareto[j][k] <= pareto[i][k], 1:n)
                    return j
                end
            end
        end
    end
    return nothing
end

function is_nondominated(
    y :: Vector{Int},
    pareto :: Vector{Vector{Int}};
    maximization :: Bool = true
) :: Bool
    if isempty(pareto)
        return true
    end
    n = length(y)
    for p in pareto
        if maximization
            if all(i -> y[i] <= p[i], 1:n)
                return false
            end
        else
            if all(i -> p[i] <= y[i], 1:n)
                return false
            end
        end
    end
    return true
end

function remove_dominated!(
    Y :: Vector{Vector{Int}},
    X :: Vector{Vector{Int}},
    y :: Vector{Int};
    maximization :: Bool = true
)
    #Remove all points of `pareto` that are dominated by `y`
    i = 1
    while i <= length(Y)
        if maximization
            if all(j -> Y[i][j] <= y[j], 1:length(y))
                deleteat!(Y, i)
                deleteat!(X, i)
            else
                i += 1
            end
        else
            if all(j -> y[j] <= Y[i][j], 1:length(y))
                deleteat!(Y, i)
                deleteat!(X, i)
            else
                i += 1
            end
        end
    end
end

function is_nondominated_set(
    pareto :: Set{Vector{Int}};
    maximization :: Bool = true
) :: Bool
    vector = pareto_vector(pareto)
    return all(
        i -> isnothing(is_nondominated(i, vector, maximization=maximization)),
        1:length(pareto)
    )
end

"""
For debugging.
Prints all dominated points in `pareto`.
"""
function print_dominated(
    pareto :: Vector{Vector{Int}};
    maximization :: Bool = true
)
    for i in 1:length(pareto)
        dominated = is_nondominated(i, pareto, maximization = maximization)
        if !isnothing(dominated)
            println("Solution ", i, " = ", pareto[i], " dominated by ", dominated, " = ", pareto[dominated])
        end
    end
end

function print_dominated(
    pareto :: Set{Vector{Int}};
    maximization :: Bool = true
)
    print_dominated(pareto_vector(pareto), maximization=maximization)
end

"""
Checks whether `pareto1` is contained by `pareto2`.
"""
function contained_by(
    pareto1 :: Vector{Vector{Int}},
    pareto2 :: Vector{Vector{Int}}
) :: Bool
    pareto2_set = Set{Vector{Int}}(pareto2)
    for nondominated in pareto1
        if !(nondominated in pareto2_set)
            return false
        end
    end
    return true
end

function contained_by(
    pareto1 :: Set{Vector{Int}},
    pareto2 :: Set{Vector{Int}}
) :: Bool
    return contained_by(
        pareto_vector(pareto1), pareto_vector(pareto2)
    )
end

"""
Changes direction of the given Pareto set by changing signs (i.e., transform
a solution from a minimization problem into maximization and vice-versa).
"""
function change_direction!(
    pareto :: Union{Vector{Vector{Int}}, Set{Vector{Int}}}
)
    for solution in pareto
        n = length(solution)
        for i in 1:n
            solution[i] = -solution[i]
        end
    end
end

"""
    is_feasible(x :: Vector{Int}, A :: Matrix{Int}, b :: Vector{Int})

    Check whether `x` is feasible for `Ax = b, x \\geq 0`.
"""
function is_feasible(x :: Vector{Int}, A :: Matrix{Int}, b :: Vector{Int}) :: Bool
    n = length(x)
    @assert size(A, 1) == length(b)
    if n != size(A, 2)
        return false
    elseif any(i -> x[i] < 0, 1:n)
        return false
    end
    return A * x == b
end

function is_efficient_set_feasible(
    efficient_set :: Vector{Vector{Int}},
    A :: Matrix{Int},
    b :: Vector{Int}
) :: Bool
    return all(x -> is_feasible(x, A, b), efficient_set)
end

function ideal_point(
    instance :: IPInstance,
    optimizer :: DataType = IPGBs.DEFAULT_SOLVER
) :: Vector{Int}
    model, vars, _ = IPInstances.jump_model(instance, optimizer=optimizer)
    sense = Int(objective_sense(model))
    ideal = Int[]
    for i in 1:num_objectives(instance)
        @objective(model, OptimizationSense(sense), instance.C[i, :]' * vars)
        optimize!(model)
        if termination_status(model) != MOI.OPTIMAL
            throw(ArgumentError("Ideal point is unbounded"))
        end
        push!(ideal, round(Int, objective_value(model)))
    end
    @objective(model, OptimizationSense(sense), instance.C[1, :]' * vars)
    return ideal
end

function nadir_bound(
    instance :: IPInstance,
    optimizer :: DataType = IPGBs.DEFAULT_SOLVER
) :: Vector{Int}
    model, vars, _ = IPInstances.jump_model(instance, optimizer=optimizer)
    sense = Int(objective_sense(model))
    nadir = Int[]
    for i in 1:num_objectives(instance)
        @objective(model, OptimizationSense(1 - sense), instance.C[i, :]' * vars)
        optimize!(model)
        if termination_status(model) != MOI.OPTIMAL
            throw(ArgumentError("Nadir bound is unbounded"))
        end
        push!(nadir, round(Int, objective_value(model)))
    end
    @objective(model, OptimizationSense(sense), instance.C[1, :]' * vars)
    return nadir
end

end
