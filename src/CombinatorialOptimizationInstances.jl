module CombinatorialOptimizationInstances

using JuMP
import LinearAlgebra
using Random

using IPGBs.GBTools
using IPGBs.IPInstances

export generate_knapsack, generate_lap, generate_set_cover, generate_set_packing
export generate_all_instances

#Knapsack problems
function generate_knapsack(
    n :: Int,
    m :: Int = 1;
    binary :: Bool = false,
    correlation :: Bool = false,
    max_coef :: Int = 1000
)
    eps = round.(Int, max_coef / 10)
    C = rand(1:max_coef, 1, n)
    A = Matrix{Int}(undef, m, n)
    for i in 1:m
        for j in 1:n
            if correlation
                #Positive correlation between values and weights makes the problem harder,
                #in average
                A[i, j] = rand((C[j] - eps):(C[j] + eps))
            else
                A[i, j] = rand(1:max_coef)
            end
        end
    end
    b = round.(Int, sum(A, dims=2)[:, 1] / 2)
    model = Model()
    if binary
        @variable(model, x[1:n] >= 0, Bin)
    else
        @variable(model, x[1:n] >= 0, Int)
    end
    @objective(model, Max, sum(C[1, j] * x[j] for j in 1:n))
    @constraint(model, [i in 1:m], sum(A[i, j] * x[j] for j in 1:n) <= b[i])
    return model, x
end

#Linear assignment problems
function generate_lap(n :: Int)
    obj = rand(1:n, n, n)
    model = Model()
    @variable(model, x[1:n, 1:n], Bin)
    @objective(model, Min, sum(obj[i, j] * x[i,j] for i in 1:n, j in 1:n))
    for j in 1:n
        @constraint(model, sum(x[i,j] for i in 1:n) == 1)
    end
    for i in 1:n
        @constraint(model, sum(x[i,j] for j in 1:n) == 1)
    end
    return model, x
end

function repeated_subsets(subsets :: Vector{Set{Int}})
    return any(subsets[i] == subsets[j] for i in 1:length(subsets), j in 1:length(subsets) if i != j)
end

function is_feasible_set_cover(subsets :: Vector{Set{Int}}, n :: Int)
    if isempty(subsets)
        return n == 0
    end
    return foldl(union, subsets) == Set(1:n)
end

is_feasible_set_packing(s :: Vector{Set{Int}}, :: Int) = length(s) > 0

function generate_subsets(n :: Int, m :: Int, p :: Float64, feasibility_check)
    subsets = Vector{Set{Int}}()
    while !feasibility_check(subsets, n) || repeated_subsets(subsets)
        empty!(subsets)
        while length(subsets) < m
            subset = Set{Int}()
            for j in 1:n
                if rand() < p
                    push!(subset, j)
                end
            end
            if isempty(subset)
                continue
            end
            push!(subsets, subset)
        end
    end
    return subsets
end

#Set covering problems
#Given n elements (say, the numbers 1:n) and m subsets of these elements, find
#a minimum number of subsets such that every element is contained in at least
#one of the chosen subsets.
#p is the probability any given element is in a subset.
function generate_set_cover(n :: Int, m :: Int, p :: Float64)
    subsets = generate_subsets(n, m, p, is_feasible_set_cover)
    model = Model()
    @variable(model, x[1:m], Bin)
    @objective(model, Min, sum(x[i] for i in 1:m))
    for i in 1:n
        @constraint(model, sum(x[j] for j in 1:m if i in subsets[j]) >= 1)
    end
    return model, x
end

#Set packing problems
#Given n elements (numbers 1:n) and m subsets of these elements, find a maximum
#number of subsets such that no element is contained in more than one of them.
#p is the probability any given element is in a subset.
function generate_set_packing(n :: Int, m :: Int, p :: Float64)
    subsets = generate_subsets(n, m, p, is_feasible_set_packing)
    model = Model()
    @variable(model, x[1:m], Bin)
    @objective(model, Max, sum(x[i] for i in 1:m))
    for i in 1:n
        @constraint(model, sum(x[j] for j in 1:m if i in subsets[j]) <= 1)
    end
    return model
end

function generate_bqp(n :: Int, range :: Int = 10)
    Q = rand(-range:range, n, n)
    #Set everything strictly below the diagonal to zero
    for i in 1:n
        for j in 1:(i-1)
            Q[i, j] = 0
        end
    end
    model = Model()
    @variable(model, x[1:n], Bin)
    @variable(model, z[i=1:n, j=(i+1):n], Bin)
    @objective(model, Min,
        sum(Q[i, j] * z[i, j] for i in 1:n, j in (i+1):n)
        + sum(Q[i, i] * x[i] for i in 1:n)
    )
    #Enforce that z[i, j] = x[i] * x[j] = x[i] AND x[j]
    for i in 1:n
        for j in (i+1):n
            @constraint(model, 2*z[i, j] <= x[i] + x[j])
            @constraint(model, z[i, j] >= x[i] + x[j] - 1)
        end
    end
    return model
end

function initial_lap(n)
    id = Matrix{Int}(LinearAlgebra.I, n, n)
    return reshape(id, n^2)
end

function initial_set_cover(m)
    #I have to consider slack variables
    return ones(Int, m)
end

function initial_set_packing(m)
    return zeros(Int, m)
end

###
### Generating instances for all problems
###

BASE_INSTANCE_DIR = "../instances/"

function all_binary_knapsacks(reps = 10)
    ns = [5, 10, 20, 30, 40, 50]
    full_path = BASE_INSTANCE_DIR * "knapsack_binary"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for correlated in [true, false]
        for n in ns
            for rep in 1:reps
                knapsack, _ = generate_knapsack(n, binary=true, correlation=correlated)
                corr = ""
                if correlated
                    corr = "_corr"
                end
                name = "knapsack_binary_" * string(n) * corr * "_" * string(rep) * ".mps"
                write_to_file(knapsack, full_path * "/" * name)
            end
        end
    end
end

function all_unbounded_knapsacks(reps = 10)
    ns = [10, 30, 50, 75, 100]
    full_path = BASE_INSTANCE_DIR * "knapsack_unbounded"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for correlated in [true, false]
        for n in ns
            for rep in 1:reps
                knapsack, _ = generate_knapsack(n, binary=false, correlation=correlated)
                corr = ""
                if correlated
                    corr = "_corr"
                end
                name = "knapsack_unbounded_" * string(n) * corr * "_" * string(rep) * ".mps"
                write_to_file(knapsack, full_path * "/" * name)
            end
        end
    end
end

function all_multiknapsacks(reps = 10)
    ns = [10, 15, 20, 25, 30]
    ms = [2, 3]
    full_path = BASE_INSTANCE_DIR * "knapsack_multidimensional"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for n in ns
        for m in ms
            for rep in 1:reps
                knapsack, _ = generate_knapsack(n, m)
                name = "knapsack_multidimensional_" * string(n) * "_" * string(m) * "_" * string(rep) * ".mps"
                write_to_file(knapsack, full_path * "/" * name)
            end
        end
    end
end

function all_binary_multiknapsacks(reps = 10)
    ns = [10, 15, 20, 25, 30, 35, 40, 45, 50]
    ms = [2, 3, 4, 5]
    full_path = BASE_INSTANCE_DIR * "knapsack_binary_multidimensional"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for n in ns
        for m in ms
            for rep in 1:reps
                knapsack, _ = generate_knapsack(n, m, binary=true)
                name = "knapsack_binary_multidimensional_" * string(n) * "_" * string(m) * "_" * string(rep) * ".mps"
                write_to_file(knapsack, full_path * "/" * name)
            end
        end
    end
end

function all_laps(reps = 10)
    ns = [5, 10, 15, 20]
    full_path = BASE_INSTANCE_DIR * "lap"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for n in ns
        for rep in 1:reps
            lap, _ = generate_lap(n)
            name = "lap_" * string(n) * "_" * string(rep) * ".mps"
            write_to_file(lap, full_path * "/" * name)
        end
    end
end

function all_set_covers(reps = 10)
    nmps = [
        (30, 15, 0.2), (30, 16, 0.2), (30, 17, 0.2), (30, 18, 0.2), (30, 19, 0.2), (30, 20, 0.2),
        (35, 15, 0.2), (35, 16, 0.2), (35, 17, 0.2), (35, 18, 0.2), (35, 19, 0.2), (35, 20, 0.2),
        (40, 15, 0.2), (40, 16, 0.2), (40, 17, 0.2), (40, 18, 0.2), (40, 19, 0.2), (40, 20, 0.2)
    ]
    full_path = BASE_INSTANCE_DIR * "set_cover"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for (n, m, p) in nmps
        for rep in 1:reps
            cover, _ = generate_set_cover(n, m, p)
            name = "set_cover_" * string(n) * "_" * string(m) * "_" * string(p) * "_" * string(rep) * ".mps"
            write_to_file(cover, BASE_INSTANCE_DIR * "set_cover" * "/" * name)
        end
    end
end

function all_set_packings(reps = 10)
    nmps = [
        (100, 25, 0.05), (100, 30, 0.05), (100, 35, 0.05), (100, 40, 0.05),
        (150, 25, 0.05), (150, 30, 0.05)
    ]
    full_path = BASE_INSTANCE_DIR * "set_packing"
    if !isdir(full_path)
        mkdir(full_path)
    end
    for (n, m, p) in nmps
        for rep in 1:reps
            packing = generate_set_packing(n, m, p)
            name = "set_packing_" * string(n) * "_" * string(m) * "_" * string(p) * "_" * string(rep) * ".mps"
            write_to_file(packing, BASE_INSTANCE_DIR * "set_packing" * "/" * name)
        end
    end
end

function generate_all_instances()
    Random.seed!(0)
    all_binary_knapsacks()
    all_unbounded_knapsacks()
    all_multiknapsacks()
    all_laps()
    all_set_covers()
    all_set_packings()
end

end
