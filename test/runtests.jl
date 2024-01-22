using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.GBAlgorithms
using MIPMatrixTools.GBTools
using Test

include("./Aqua.jl")
include("./test_functions.jl")

@testset "Buchberger: Binary Knapsacks" begin
    for n in [5, 10, 15, 20, 25, 30]
        println("Buchberger test for binary knapsack, n = ", n)
        ipgbs_result, fourti2_result, instance = test_buchberger(n, quiet=true)
        alg = BuchbergerAlgorithm(ipgbs_result, instance)
        truncate(b) = GBAlgorithms.general_truncate(alg, b)
        @test is_truncated_groebner_basis(current_basis(alg), truncate)
        @test length(ipgbs_result) == length(fourti2_result)
        @test GBTools.isequal(ipgbs_result, fourti2_result)
        if !GBTools.isequal(ipgbs_result, fourti2_result)
            println("4ti2 result:")
            println(fourti2_result)
            println("My result:")
            println(ipgbs_result)
        end
        println()
    end
end
