using IPGBs
using IPGBs.GBTools
using IPGBs.GBElements
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.Buchberger
using IPGBs.FourTi2
using IPGBs.SignatureAlgorithms
using IPGBs.BinomialSets
using IPGBs.IPInstances
using IPGBs.Markov

using MultiObjectiveInstances
import Random

"""
Returns the GB given by my implementation of Buchberger's algorithm and the
equivalent 4ti2 GB.
"""
function test_buchberger(
    n::Int;
    seed = 0,
    setseed = true,
    implicit_representation = false,
    truncation_type = :Heuristic,
    quiet = false
)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary = true)
    ideal_gens = [
        lattice_generator_graded(
            i, instance.A, instance.b, Float64.(instance.C), instance.u,
            check_truncation = false)
        for i in 1:size(instance.A, 2)]
    gens_4ti2 = GradedBinomials.fourti2_form(ideal_gens)

    #Use 4ti2 to compute a reduced GB for reference
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    truncate = truncation_type != :None
    trunc_sol = truncate ? initial_solution : Int[]
    rgb, time, _, _, _ = @timed groebner(
        gens_4ti2, instance_4ti2.C, truncation_sol = trunc_sol,
        lattice = true, quiet = quiet
    )
    if !quiet
        println()
        println("4ti2 results")
        @show size(rgb, 1) time
        println()
    end
    fourti2gb = IPGBs.GBTools.tovector(rgb)

    #Compute a GB using my Buchberger implementation
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, u=instance.u, use_signatures = false,
        implicit_representation = implicit_representation,
        truncation_type = truncation_type, quiet = quiet
    )
    if !quiet
        println()
        println("my results")
        @show length(gb) time
    end
    return gb, fourti2gb
end

function test_siggb(
    n::Int;
    seed = 0,
    setseed = true,
    module_order = :ltpot,
    truncation_type = :Heuristic,
    quiet = false
)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}},Vector{Vector{Int}}}
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary = true)
    lattice_basis = [
        lattice_generator_graded(
            i, instance.A, instance.b, Float64.(instance.C), instance.u,
            check_truncation = false)
        for i in 1:size(instance.A, 2)]
    lattice_4ti2 = GradedBinomials.fourti2_form(lattice_basis)

    #Use 4ti2 to compute a reduced GB for reference
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
        instance_4ti2)
    truncate = truncation_type != :None
    trunc_sol = truncate ? initial_solution : Int[]
    rgb, time, _, _, _ = @timed groebner(
        lattice_4ti2, instance_4ti2.C, truncation_sol = trunc_sol,
        lattice = true, quiet = quiet
    )
    if !quiet
        println()
        println("4ti2 results")
        @show size(rgb, 1) time
        println()
    end
    fourti2gb = IPGBs.GBTools.tovector(rgb)

    #My results
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u, use_signatures = true,
        module_order = module_order, truncation_type = truncation_type, quiet = true
    )
    if !quiet
        println()
        println("Signature results")
        @show length(gb) time
        println()
    end
    #Basic Buchberger results
    bgb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, instance.u, use_signatures = false,
        truncation_type = truncation_type, quiet = true
    )
    if !quiet
        println()
        println("Buchberger results")
        @show length(bgb) time
    end
    return gb, fourti2gb, bgb
end

function run_algorithm(
    n::Int;
    seed = 0,
    setseed = true,
    module_order = :ltpot,
    use_signatures = false,
    implicit_representation = false,
    truncation_type = :Heuristic,
    quiet = false,
    minimization = true,
    binary=true
)
    if setseed
        Random.seed!(seed)
    end
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary = binary)
    #My results
    gb, time, _, _, _ = @timed groebner_basis(
        instance.A, instance.b, instance.C, u=instance.u,
        use_signatures = use_signatures, module_order = module_order,
        implicit_representation = implicit_representation,
        truncation_type = truncation_type, quiet = quiet,
        minimization = minimization
    )
    if !quiet
        println()
        println("My results")
        @show length(gb) time
    end
    return gb
end

function test_markov(
    n::Int;
    seed = 0,
    setseed = true
)::Vector{Vector{Int}}
    if setseed
        Random.seed!(seed)
    end
    #TODO should probably generate harder instances for Markov bases, as these
    #knapsacks lead to very small bases...
    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary = true)
    #Compute Markov basis with my implementation of project-and-lift
    mb = markov_basis(instance.A, instance.b, instance.C, instance.u,
        algorithm = :ProjectAndLift)
    #Compute (truncated) Markov basis with 4ti2's implementation
    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
    initsol = MultiObjectiveInstances.Knapsack.knapsack_initial(instance_4ti2)
    correct_mb = markov(instance_4ti2.A, instance_4ti2.C, truncation_sol = initsol)
    @show mb correct_mb
    return mb
end
