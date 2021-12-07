"""
Experiment with cover inequalities of binary knapsacks. The idea is checking what
happens to GBs (and the performance of Buchberger's algorithm) as I add cover
inequalities.

Initial results:
- the GBs do get bigger with the cover inequalities
- the running time also gets a lot larger... probably due to the GB sizes
- my main hypothesis for this is that the cover inequalities make the < b fibers
more complicated, leading to more stuff in the GB
- adding a single cover inequality may or may not change the GB size
- I have no idea why some cover inequalities make the basis larger and others
don't. I can't see a pattern, either, checking test_valid(5)

TODO move these cover inequality generating functions to
MultiObjectiveInstances.Knapsacks
"""

using IPGBs
using IPGBs.FourTi2

using MultiObjectiveInstances
import Combinatorics
import Random

"""
Generates all cover inequalities for a binary knapsack.
Each is represented by a vector of variable indices occurring in the
inequality.

If max_inequalities != nothing, generates only max_inequalities cover
inequalities instead.
"""
function cover_inequalities(instance, n; max_inequalities=nothing)
    a = instance.A[1, :]
    ineqs = []
    if !isnothing(max_inequalities) && max_inequalities == 0
        return ineqs
    end
    for var_subset in Combinatorics.powerset(1:n)
        if !isempty(var_subset)
            s = sum(a[i] for i in var_subset)
            if s > instance.b[1]
                push!(ineqs, var_subset)
                if !isnothing(max_inequalities) && length(ineqs) == max_inequalities
                    return ineqs
                end
            end
        end
    end
    return ineqs
end

"""
Generates instance data for `instance` that includes cover inequalities.
"""
function add_cover_inequalities(instance, n; max_inequalities=nothing)
    ineqs = cover_inequalities(instance, n, max_inequalities=max_inequalities)
    new_vars = length(ineqs) #One new slack per inequality
    total_vars = size(instance.A, 2) + new_vars
    new_C = [instance.C zeros(Int, size(instance.C, 1), new_vars)]
    #Update A, b
    new_b = copy(instance.b)
    new_A = [instance.A zeros(Int, size(instance.A, 1), new_vars)]
    j = 1
    for ineq in ineqs
        new_constr = zeros(Int, 1, total_vars)
        for index in ineq
            new_constr[1, index] = 1
        end
        new_constr[1, size(instance.A, 2) + j] = 1
        push!(new_b, length(ineq) - 1)
        new_A = [new_A; new_constr]
        j += 1
    end
    initial_sol = [zeros(Int, n); new_b ]
    return new_A, new_b, new_C, initial_sol
end

function test_cover_one_by_one(instance, n; max_inequalities=nothing)
    @show instance
    println("Trying cover inequalities one by one")
    ineqs = cover_inequalities(instance, n, max_inequalities=max_inequalities)
    new_vars = 1
    total_vars = size(instance.A, 2) + new_vars
    new_C = [instance.C zeros(Int, size(instance.C, 1), new_vars)]
    new_b = [instance.b; 0]
    new_A = [instance.A zeros(Int, size(instance.A, 1), new_vars); zeros(Int, 1, total_vars)]
    for ineq in ineqs
        #Compute new matrix and RHS including the cover ineq
        new_b[length(new_b)] = length(ineq) - 1
        for index in 1:size(new_A, 2)
            new_A[size(new_A, 1), index] = 0
        end
        for index in ineq
            new_A[size(new_A, 1), index] = 1
        end
        new_A[size(new_A, 1), size(new_A, 2)] = 1
        initial_sol = [ zeros(Int, n); new_b ]

        #Compute the GB with the new cover ineq
        @show ineq new_b[length(new_b)]
        gb, t, _, _, _ = @timed groebner(new_A, new_C, truncation_sol=initial_sol)
        @show size(gb, 1) t
    end
end

function test_valid(
    n :: Int,
    seed = 0,
    setseed = true;
    max_inequalities = nothing
)
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true, format="4ti2")
    initial_sol = MultiObjectiveInstances.Knapsack.knapsack_initial(instance)
    #Compute the original GB
    gb, t, _, _, _ = @timed groebner(instance.A, instance.C, truncation_sol=initial_sol)
    @show size(gb, 1) t

    #Compute the GB with cover inequalities
    A, b, C, initial_sol2 = add_cover_inequalities(instance, n, max_inequalities=max_inequalities)
    gb2, t2, _, _, _ = @timed groebner(A, C, truncation_sol=initial_sol2)
    @show size(gb2, 1) t2

    test_cover_one_by_one(instance, n, max_inequalities=max_inequalities)
end
