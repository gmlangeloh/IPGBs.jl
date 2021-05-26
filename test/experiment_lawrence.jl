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

Simple results:
- the theorem holds in practice (obviously!). Now I can better visualize it
- I didn't compute the GBs wrt many objective functions, it is unnecessary. It will
work, when the "full" Lawrence lifting is done
- dropping the last lifting constraint (corresponding to the slack variable of the
knapsack) makes a huge difference. The GB is already nowhere near a Graver basis
when this constraint is removed
- however, dropping more constraints keeps making the problem easier. So, if there
was a way to deal with them algorithmically, it might be worth it regardless.
- this trend continues even if we are using truncated GBs instead of full GBs.
So, an algorithm that deals with upper bound constraints separately sounds useful
even for this case (which is what I expected, but it is good to be sure...)
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

function partial_lifting(
    A :: Array{Int, 2},
    n_lift :: Int
) :: Array{Int, 2}
    m, n = size(A)
    Zmnlift = zeros(Int, m, n_lift)
    Inlift = I(n_lift)
    Znlift = zeros(Int, n_lift, n - n_lift)
    return [A Zmnlift; Inlift Znlift Inlift]
end

function test_lawrence(
    n :: Int,
    seed = 0,
    setseed = true;
    check_graver = false,
    max_lifting = n,
    truncate = false
)
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, format="4ti2")
    #To truncate, I need to compute a feasible solution of each problem
    sol = []
    if truncate
        sol = [ zeros(Int, n); instance.b ]
    end

    println("Before Lawrence lifting:")
    gb = groebner(instance.A, truncation_sol=sol)
    println("GB size: ", size(gb, 1))
    if check_graver
        gra = graver(instance.A)
        println("Graver size: ", size(gra, 1))
    end

    println("After Lawrence lifting")
    for n_lift in max_lifting:-1:1
        new_sol = []
        if truncate
            if n_lift == size(instance.A, 2)
                new_sol = [ sol; ones(Int, n_lift); 0]
            else
                new_sol = [ sol; ones(Int, n_lift) ]
            end
        end
        println("Adding bounds for first ", n_lift, " variables")
        partial_A = partial_lifting(instance.A, n_lift)
        gb = groebner(partial_A, truncation_sol=new_sol)
        println("GB size: ", size(gb, 1))
        if check_graver
            gra = graver(partial_A)
            println("Graver size: ", size(gra, 1))
        end
    end
end
