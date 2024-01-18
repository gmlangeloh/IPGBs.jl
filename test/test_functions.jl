using IPGBs
using IPGBs.FourTi2

using MIPMatrixTools.GBTools
using MIPMatrixTools.CombinatorialOptimizationInstances

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
    instance, _ = generate_knapsack(n, binary = true)
    rgb, time, _, _, _ = @timed groebner(instance)
    if !quiet
        println()
        println("4ti2 results")
        @show size(rgb, 1) time
        println()
    end
    fourti2gb = IPGBs.GBTools.tovector(rgb)

    #Compute a GB using my Buchberger implementation
    gb, time, _, _, _ = @timed groebner_basis(
        instance, use_signatures = false,
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

#TODO: Reintroduce this function after updating the signature algorithm
#function test_siggb(
#    n::Int;
#    seed = 0,
#    setseed = true,
#    module_order = :ltpot,
#    truncation_type = :Heuristic,
#    quiet = false
#)::Tuple{Vector{Vector{Int}},Vector{Vector{Int}},Vector{Vector{Int}}}
#    if setseed
#        Random.seed!(seed)
#    end
#    instance = MultiObjectiveInstances.Knapsack.knapsack_A(n, binary = true)
#    lattice_basis = [
#        lattice_generator_graded(
#            i, instance.A, instance.b, Float64.(instance.C), instance.u,
#            check_truncation = false)
#        for i in 1:size(instance.A, 2)]
#    lattice_4ti2 = GradedBinomials.fourti2_form(lattice_basis)
#
#    #Use 4ti2 to compute a reduced GB for reference
#    instance_4ti2 = MultiObjectiveInstances.fourti2_stdform(instance)
#    initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(
#        instance_4ti2)
#    truncate = truncation_type != :None
#    trunc_sol = truncate ? initial_solution : Int[]
#    rgb, time, _, _, _ = @timed groebner(
#        lattice_4ti2, instance_4ti2.C, truncation_sol = trunc_sol,
#        lattice = true, quiet = quiet
#    )
#    if !quiet
#        println()
#        println("4ti2 results")
#        @show size(rgb, 1) time
#        println()
#    end
#    fourti2gb = IPGBs.GBTools.tovector(rgb)
#
#    #My results
#    gb, time, _, _, _ = @timed groebner_basis(
#        instance.A, instance.b, instance.C, instance.u, use_signatures = true,
#        module_order = module_order, truncation_type = truncation_type, quiet = true
#    )
#    if !quiet
#        println()
#        println("Signature results")
#        @show length(gb) time
#        println()
#    end
#    #Basic Buchberger results
#    bgb, time, _, _, _ = @timed groebner_basis(
#        instance.A, instance.b, instance.C, instance.u, use_signatures = false,
#        truncation_type = truncation_type, quiet = true
#    )
#    if !quiet
#        println()
#        println("Buchberger results")
#        @show length(bgb) time
#    end
#    return gb, fourti2gb, bgb
#end
