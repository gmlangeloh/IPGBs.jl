using IPGBs
using IPGBs.IPInstances
using IPGBs.Markov
import IPGBs.FourTi2

using JuMP
using Random
using GLPK

function generate_random_instance(n :: Int, m :: Int, coef_range :: Int)
    Random.seed!(0)
    feasible = false
    instance = nothing
    while !feasible
        A = rand(-coef_range:coef_range, m, n)
        b = rand((-2*coef_range):(2*coef_range), m)
        C = rand(-coef_range:coef_range, n)

        model = Model(GLPK.Optimizer)
        @variable(model, x[1:n] >= 0, Int)
        @constraint(model, A*x .<= b)
        @objective(model, Max, C' * x)
        instance = IPInstance(model)
        optimize!(model)
        feasible = termination_status(model) == MOI.OPTIMAL 
    end
    return instance
end

function run_random_instance(n :: Int, m :: Int, coef_range :: Int = 10)
    instance = generate_random_instance(n, m, coef_range)
    println("Solving instance with n = $n, m = $m, coef_range = $coef_range")
    println("Instance: ", instance)
    println("Solving with IPGBs")
    println("First, I will compute a Markov basis")
    markov = markov_basis(instance)
    println("Markov basis size:", length(markov))
    println(markov)
    gb = groebner_basis(instance)
    println("GB size:", length(gb))
    println(gb)
    println("Solving with 4ti2")
    println("First, I will compute a Markov basis")
    markov2 = FourTi2.markov(instance)
    println("Markov basis size:", size(markov2, 1))
    println(markov2)
    println("Now, I will compute a Gröbner basis")
    gb2 = FourTi2.groebner(instance)
    println("GB size:", size(gb2, 1))
    println(gb2)
    println("Is equal? ", length(gb) == size(gb2, 1))
end