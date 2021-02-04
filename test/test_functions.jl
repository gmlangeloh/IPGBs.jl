using IPGBs
using IPGBs.GBElements
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.Buchberger
using IPGBs.FourTi2
using IPGBs.SignatureAlgorithms

using MultiObjectiveInstances
import Random

"""
Returns the GB given by my implementation of Buchberger's algorithm and the
equivalent 4ti2 GB.
"""
function test_buchberger(
    n :: Int;
    seed = 0,
    setseed = true,
    implicit_representation = false,
    truncate = true
) :: Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true)
    lattice_basis = [
        lattice_generator_graded(i, instance.A, instance.b, instance.C, instance.u,
                                 check_truncation=false)
        for i in 1:size(instance.A, 2)]
    lattice_4ti2 = fourti2_form(lattice_basis)

    #Use 4ti2 to compute a reduced GB for reference
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    trunc_sol = truncate ? initial_solution : Int[]
    rgb, time, _, _, _ = @timed groebner(
        lattice_4ti2, instance_4ti2.C, truncation_sol=trunc_sol,
        lattice=true, quiet=false
    )
    println()
    println("4ti2 results")
    @show size(rgb, 1) time
    fourti2gb = IPGBs.GBTools.tovector(rgb)
    println()

    #Compute a GB using my Buchberger implementation
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u, use_signatures=false,
        implicit_representation=implicit_representation,
        truncate=truncate
    )
    println()
    println("my results")
    @show length(gb) time

    return gb, fourti2gb
end

function test_siggb(
    n :: Int;
    seed = 0,
    setseed = true,
    module_order = :ltpot,
    truncate = true
) :: Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Vector{Int}}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true)
    lattice_basis = [
        lattice_generator_graded(i, instance.A, instance.b, instance.C, instance.u,
                                 check_truncation=false)
        for i in 1:size(instance.A, 2)]
    lattice_4ti2 = fourti2_form(lattice_basis)

    #Use 4ti2 to compute a reduced GB for reference
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    trunc_sol = truncate ? initial_solution : Int[]
    rgb, time, _, _, _ = @timed groebner(
        lattice_4ti2, instance_4ti2.C, truncation_sol=trunc_sol,
        lattice=true
    )
    println()
    println("4ti2 results")
    @show size(rgb, 1) time
    fourti2gb = IPGBs.GBTools.tovector(rgb)
    println()

    #My results
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u, use_signatures=true,
        module_order=module_order, truncate=truncate
    )
    println()
    println("Signature results")
    @show length(gb) time
    println()
    #Basic Buchberger results
    bgb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u, use_signatures=false,
        truncate=truncate
    )
    println()
    println("Buchberger results")
    @show length(bgb) time
    return gb, fourti2gb, bgb
end

function run_algorithm(
    n :: Int;
    seed = 0,
    setseed = true,
    module_order = :ltpot,
    use_signatures = true,
    implicit_representation = false,
    truncate = true
) :: Vector{Vector{Int}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary=true)
    lattice_basis = [
        lattice_generator_graded(i, instance.A, instance.b, instance.C, instance.u,
                                 check_truncation=false)
        for i in 1:size(instance.A, 2)]
    #My results
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u,
        use_signatures=use_signatures, module_order=module_order,
        implicit_representation=implicit_representation,
        truncate=truncate
    )
    println()
    println("My results")
    @show length(gb) time
    return gb
end
