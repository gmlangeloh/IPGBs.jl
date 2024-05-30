module Walkback

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Optimize
using IPGBs.SupportTrees
using MIPMatrixTools.IPInstances
using MIPMatrixTools.SolverTools

using JuMP

export enumerate_solutions

function enumerate_solutions(
    instance :: IPInstance
) :: Set{Vector{Int}}
    gb = groebner_basis(instance)
    return enumerate_solutions(instance, gb)
end

function enumerate_solutions(
    instance :: IPInstance,
    gb :: Vector{Vector{Int}}
) :: Set{Vector{Int}}
    binomial_gb = BinomialSet(gb, instance, Binomial)
    return enumerate_solutions(instance, binomial_gb)
end

function enumerate_solutions(
    instance :: IPInstance,
    gb :: BinomialSet
) :: Set{Vector{Int}}

    #Step 0: Verify that the polyhedron of `instance` is bounded
    if !SolverTools.is_bounded_polyhedron(instance.A)
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
    feasible_solutions = Set{Vector{Int}}()
    solution_queue = Vector{Vector{Int}}()
    push!(solution_queue, opt_sol)
    while !isempty(solution_queue)
        current_sol = popfirst!(solution_queue)
        if current_sol in feasible_solutions
            continue
        end
        push!(feasible_solutions, current_sol)
        reducers = enumerate_reducers(current_sol, gb, reduction_tree(gb), negative=true)
        for reducer in reducers
            bin_sol = to_gbelement(current_sol, order(gb), Binomial, false)
            GBElements.reduce!(bin_sol, reducer, negative=true)
            new_sol = copy(element(bin_sol))
            if !(new_sol in feasible_solutions)
                push!(solution_queue, new_sol)
            end
        end
    end
    return feasible_solutions
end

end
