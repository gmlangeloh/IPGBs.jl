using IPGBs
using IPGBs.CombinatorialOptimizationInstances
using IPGBs.FourTi2
using IPGBs.GBTools
using IPGBs.IPInstances

using JuMP
import Random

function compare_to_solver(ipgbs_value, ipgbs_status, solver_value, solver_status)
    #If both status are OPTIMAL, compare the values
    #Otherwise, pass if both status are the same
    same_status = ipgbs_status == solver_status
    if !same_status
        return false
    end
    if ipgbs_status == MOI.OPTIMAL
        return ipgbs_value == solver_value
    end
    return true #Other status should pass automatically
end

"""
Returns the GB given by my implementation of Buchberger's algorithm and the
equivalent 4ti2 GB.
"""
function ipgbs_and_fourti2(
    n::Int;
    seed = 0,
    setseed = true,
    implicit_representation = false,
    truncation_type = :Heuristic,
    quiet = false
)
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
    fourti2gb = GBTools.tovector(rgb)

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
    return gb, fourti2gb, IPInstance(instance)
end

function test_buchberger_correctness(ipgbs_result, fourti2_result, instance)
    alg = BuchbergerAlgorithm(ipgbs_result, instance)
    truncate(b) = GBAlgorithms.general_truncate(alg, b)
    println("Is truncated GB? ", @test is_truncated_groebner_basis(current_basis(alg), truncate))
    println("4ti2 and IPGBs have GBs of same length? ", @test length(ipgbs_result) == length(fourti2_result))
    println("4ti2 and IPGBs return the same GB? ", @test GBTools.isequal(ipgbs_result, fourti2_result))
    println()
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
