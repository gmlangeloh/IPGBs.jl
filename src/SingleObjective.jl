"""
Tools for solving single objective integer programs given by their model
min c' * x
s.t. A * x = b
x >= 0
"""
module SingleObjective
export usesolver, linear_relaxation, feasible_solution

using MultiObjectiveInstances

using JuMP
using GLPK

"""
Creates a JuMP model (for the given solver) from `A`, `b` and `c`.
"""
function makemodel(
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Vector{Int};
    optimizer = GLPK.Optimizer,
    integer :: Bool = true
) :: Tuple{JuMP.Model, Array{JuMP.VariableRef}}
    model = Model(optimizer)
    set_silent(model)
    n = length(c)
    if integer
        @variable(model, x[1:n] >= 0, Int)
    else
        @variable(model, x[1:n] >= 0)
    end
    obj = c' * x
    @objective(model, Min, obj)
    for i in 1:size(A, 1)
        row = A[i, :]
        @constraint(model, row' * x == b[i])
    end
    return model, x
end

"""
Uses JuMP / Cbc to solve an integer programming problem given in the form
min c' * x
s.t. A * x = b
x >= 0
"""
function usesolver(
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Vector{Int};
    initial_solution = Int[],
    optimizer = GLPK.Optimizer
) :: Tuple{Vector{Int}, Int}
    model, x = makemodel(A, b, c, optimizer=optimizer)
    if length(initial_solution) > 0
        for i in eachindex(initial_solution)
            set_start_value(x[i], initial_solution[i])
        end
    end
    optimize!(model)
    solution = value.(x)
    objective = objective_value(model)
    return Int.(round.(solution)), Int(round(objective))
end

"""
    lexmin(
    A :: Matrix{Int},
    b :: Vector{Int},
    C :: Matrix{Int},
    i :: Int;
    solver :: String = "Cbc"
) :: Tuple{Vector{Int}, Int}

Solve the lexicographical minimization of the multiobjective problem

min. C * x
s.t. A * x == b
x >= 0, x integer

considering the i-th objective as the first objective, and optimizing
the remaining objectives as tiebreakers in order 1, 2, ..., i - 1,
i + 1, ... p.
"""
function lexmin(
    A :: Matrix{Int},
    b :: Vector{Int},
    C :: Matrix{Int},
    i :: Int;
    optimizer = GLPK.Optimizer
) :: Vector{Int}
    objective = C[i, :]
    model, x = makemodel(A, b, objective, optimizer=optimizer)
    optimize!(model)
    val = round(Int, objective_value(model))
    for j in 1:size(C, 1)
        if j != i
            @constraint(model, objective' * x == val)
            objective = C[j, :]
            @objective(model, Min, objective)
            optimize!(model)
            val = round(Int, objective_value(model))
        end
    end
    return round.(Int, value.(x))
end

function feasible_solution(
    A :: Matrix{Int},
    b :: Vector{Int};
    optimizer = GLPK.Optimizer
) :: Vector{Int}
    n = size(A, 2)
    c = zeros(Int, n)
    solution, _ = usesolver(A, b, c, optimizer=optimizer)
    return solution
end

"""
Solves the linear relaxation of the integer program
min c' * x
s.t. A * x = b
x >= 0
using JuMP.
"""
function linear_relaxation(
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Vector{Int};
    optimizer = GLPK.Optimizer
) :: Tuple{Vector{Float64}, Float64}
    model, x = makemodel(A, b, c, optimizer=optimizer, integer=false)
    optimize!(model)
    solution = value.(x)
    objective = objective_value(model)
    return solution, objective
end

end
