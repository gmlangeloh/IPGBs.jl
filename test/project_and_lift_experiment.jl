using MultiObjectiveInstances.Knapsack
using IPGBs
using IPGBs.IPInstances
using IPGBs.Markov

knapsack = knapsack_A(5, format="4ti2", binary=true)
u = ones(Int, size(knapsack.A, 2))
instance = IPInstances.IPInstance(
    knapsack.A, knapsack.b, knapsack.C, u, apply_normalization=false,
    invert_objective=false)
mbasis = Markov.project_and_lift(instance)
@show mbasis
