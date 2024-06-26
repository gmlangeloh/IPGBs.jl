module Walkback

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.IPInstances
using IPGBs.Optimize
using IPGBs.SolverTools
using IPGBs.SupportTrees

using JuMP

export enumerate_solutions

function enumerate_solutions(
    instance :: IPInstance;
    max_solutions :: Int = 0
) :: Set{Vector{Int}}
    gb = groebner_basis(instance)
    return enumerate_solutions(instance, gb, max_solutions=max_solutions)
end

function enumerate_solutions(
    instance :: IPInstance,
    gb :: Vector{Vector{Int}};
    max_solutions :: Int = 0
) :: Set{Vector{Int}}
    binomial_gb = BinomialSet(gb, instance, Binomial)
    return enumerate_solutions(instance, binomial_gb, max_solutions=max_solutions)
end

function enumerate_solutions(
    instance :: IPInstance,
    gb :: BinomialSet; #Inverted GB
    max_solutions :: Int = 0
) :: Set{Vector{Int}}

    #Step 0: Verify that the polyhedron of `instance` is bounded
    if !SolverTools.is_bounded_polyhedron(instance.A, optimizer=instance.optimizer)
        throw(ArgumentError(
            "The polyhedron of `instance` is unbounded,
            there may be infinite feasible solutions.")
        )
    end

    #Step 1: get optimal solution of `instance`
    opt_sol, _, status = IPInstances.solve(instance)
    if status != MOI.OPTIMAL
        return Set{Vector{Int}}()
    end
    optimize_with!(opt_sol, gb) #Tiebreak using `gb`

    #Step 2: Walk backwards through the gb-skeleton
    inv_gb = BinomialSets.inverse_set(gb)
    feasible_solutions = Set{Vector{Int}}()
    solution_queue = Vector{Vector{Int}}()
    push!(solution_queue, opt_sol)
    while !isempty(solution_queue) &&
        (max_solutions == 0 || length(feasible_solutions) < max_solutions)

        current_sol = popfirst!(solution_queue)
        if current_sol in feasible_solutions
            continue
        end
        push!(feasible_solutions, current_sol)
        reducers = enumerate_reducers(current_sol, inv_gb, reduction_tree(inv_gb))
        for reducer in reducers
            bin_sol = to_gbelement(current_sol, order(gb), Binomial, false)
            result = bin_sol - reducer
            new_sol = copy(element(result))
            if !(new_sol in feasible_solutions)
                push!(solution_queue, new_sol)
            end
        end
    end
    return feasible_solutions
end

end
