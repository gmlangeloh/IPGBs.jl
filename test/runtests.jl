using IPGBs
using IPGBs.GBElements
using IPGBs.GradedBinomials
using IPGBs.Buchberger
using IPGBs.FourTi2
using Test
using MultiObjectiveInstances

import Random

function test_buchberger(
    n :: Int;
    seed = 0,
    setseed = true,
    structure = Binomial
) :: Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true)
    lattice_basis = [ lattice_generator_graded(i, instance.A, instance.C)
                      for i in 1:size(instance.A, 2)]
    lattice_4ti2 = fourti2_form(lattice_basis)

    #Use 4ti2 to compute a reduced GB for reference
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    rgb, time, _, _, _ = @timed groebner(
        lattice_4ti2, instance_4ti2.C, truncation_sol=initial_solution,
        lattice=true
    )
    println("4ti2 results")
    @show size(rgb, 1) time
    fourti2gb = IPGBs.GBTools.tovector(rgb)

    #Compute a GB using my Buchberger implementation
    gb, time, _, _, _ = @timed buchberger(
        instance.A, instance.b, instance.C, instance.u, structure=structure
    )
    println("my results")
    @show length(gb) time

    return gb, fourti2gb
end

function test_siggb(
    n :: Int;
    seed = 0,
    setseed = true
)
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true)
    #4ti2 stuff
    instance_4ti2 = fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    rgb, time, _, _, _ = @timed FourTi2.groebner(
        instance_4ti2.A, instance_4ti2.C, truncation_sol=initial_solution
    )
    @show instance.C
    println("4ti2 results")
    @show size(rgb, 1) time
    #My results
    gb, time, _, _, _ = @timed siggb(instance.A, instance.b, instance.C,
                                     instance.u)
    println("Signature results")
    @show length(gb) time
    #Basic Buchberger results
    bgb, time, _, _, _ = @timed buchberger(instance)
    println("Buchberger results")
    @show length(bgb) time
    return gb, rgb, bgb
end

@testset "IPGBs.jl" begin
    for s in [ Binomial, GradedBinomial ]
        for n in [5, 10, 15, 20, 25]
            println("Buchberger test for ", s, " structure, n = ", n)
            gb, fourti2gb = test_buchberger(n, structure=s)
            @test IPGBs.GBTools.isequal(gb, fourti2gb)
            println()
        end
    end

    #TODO Add tests for siggb
end
