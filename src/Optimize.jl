module Optimize

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.GBElements

using MIPMatrixTools.IPInstances

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

end
