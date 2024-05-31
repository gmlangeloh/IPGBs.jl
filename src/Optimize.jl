module Optimize

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Markov

using MIPMatrixTools.IPInstances

using JuMP

export optimize_with!

"""
    optimize_with(
    solution :: Vector{Int},
    instance :: IPInstance,
    test_set :: Vector{Vector{Int}}
) :: Vector{Int}

    Use `test_set` to optimize `solution`, a feasible solution to `instance`,
    returning a new feasible solution that cannot be further improved by
    `test_set`.

    If `test_set` is a Gröbner basis, then the returned solution is optimal.
    Otherwise, it is simply a locally optimal solution with respect to `test_set`.
"""
function optimize_with!(
    solution :: Vector{Int},
    instance :: IPInstance,
    test_set :: Vector{Vector{Int}}
)
    binomial_gb = BinomialSet(test_set, instance, Binomial)
    optimize_with!(solution, instance, binomial_gb)
end

function optimize_with!(
    solution :: Vector{Int},
    test_set :: BinomialSet
)
    bin_sol = to_gbelement(solution, order(test_set), Binomial, false)
    BinomialSets.reduce!(bin_sol, test_set, is_monomial_reduction = true)
    copyto!(solution, element(bin_sol))
end

function optimize(
    instance :: IPInstance;
    completion :: Symbol = :IPGBs,
    truncation_type :: Symbol = :LP,
    solution :: Vector{Int} = zeros(Int, instance.n),
    quiet :: Bool = true
) :: Tuple{Vector{Int}, Int, TerminationStatusCode}
    initial_solution = copy(solution)
    if !is_bounded(instance)
        return initial_solution, 0, DUAL_INFEASIBLE
    end
    if !quiet
        lp_val = IPInstances.linear_relaxation(instance)
        println("LP => ", lp_val)
    end
    _, is_opt, opt_sol, opt_val = project_and_lift(
        instance, completion=completion, truncation_type=truncation_type,
        optimize=true, solution=initial_solution, quiet = quiet
    )
    @assert is_opt
    return opt_sol, opt_val, OPTIMAL
end

function optimize(
    filepath :: String;
    completion :: Symbol = :IPGBs,
    truncation_type :: Symbol = :LP,
    solution :: Vector{Int} = Int[],
    quiet :: Bool = true
) :: Tuple{Vector{Int}, Int, TerminationStatusCode}
    initial_solution = copy(solution)
    instance = IPInstance(filepath)
    if !isempty(initial_solution) && !is_feasible_solution(IPInstance(filepath), initial_solution)
        throw(ArgumentError("Initial solution is not feasible."))
    elseif isempty(initial_solution)
        initial_solution = zeros(Int, instance.n)
    end
    return optimize(
        instance, completion=completion, truncation_type=truncation_type,
        solution=initial_solution, quiet = quiet
    )
end

end
