using MultiObjectiveInstances.Knapsack
using IPGBs
using IPGBs.IPInstances
using IPGBs.Markov

knapsack = knapsack_A(5, format="4ti2", binary=true)
u = Union{Int, Nothing}[]
for i in 1:size(knapsack.A, 2)
    if i == 6 #Slack variable has no upper bound
        push!(u, nothing)
    else
        push!(u, 1)
    end
end
instance = IPInstances.IPInstance(
    knapsack.A, knapsack.b, knapsack.C, u, apply_normalization=false,
    invert_objective=false)
mbasis = Markov.project_and_lift(instance)
@show mbasis
