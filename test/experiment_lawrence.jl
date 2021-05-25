"""
Experiment with the Lawrence lifting theorem:
Let A' = [A 0; In In] be a matrix where In is the identity matrix of size n.
Then any reduced GB of I_{A'} is equal to its Universal GB and Graver basis.

In particular, I want to:
- see this happening for knapsacks, when I add the bound constraints on the variables
and the slack variable
- compute GBs wrt many objective functions
- understand what happens to these GBs when I relax each variable bound constraint

I hope this will help me understand how necessary this kind of constraint is in
practice, and whether I can simply drop them and simulate them inside the algorithm
somehow.
"""

using IPGBs
using IPGBs.FourTi2

using LinearAlgebra
using MultiObjectiveInstances
import Random

"""
Given a matrix A, computes the Lawrence lifting of A, that is,
Lambda(A) = [A 0; In In]
"""
function lawrence_lifting(
    A :: Array{Int, 2}
) :: Array{Int, 2}
    m, n = size(A)
    Zmn = zeros(Int, m, n)
    return [A Zmn; I(n) I(n)]
end

function test_lawrence(
    n :: Int,
    seed = 0,
    setseed = true
)
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, format="4ti2")
    lifted_A = lawrence_lifting(instance.A)
    #TODO lift C and find an initial solution
    #To begin with, I don't think I need to provide any of that!
    gb = groebner(lifted_A)
    println("GB size: ", size(gb, 1))
    gra = graver(lifted_A)
    println("Graver size: ", size(gra, 1))
end
