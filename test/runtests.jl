using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.FourTi2
using IPGBs.GBAlgorithms
using MIPMatrixTools.GBTools
using Test

include("./test_functions.jl")

@testset "All tests" begin

    include("./Aqua.jl")
    @testset "Buchberger: Simple binary knapsacks" begin
        for n in [5, 10, 15, 20, 25, 30]
            println("Buchberger test for binary knapsack, n = ", n)
            ipgbs_result, fourti2_result, instance = ipgbs_and_fourti2(n, quiet=true)
            test_buchberger_correctness(ipgbs_result, fourti2_result, instance)
        end
    end

    @testset "Buchberger: combinatorial optimization instances" begin
        for filename in readdir("test_instances", join=true)
            if endswith(filename, ".mps")
                println("Buchberger test for ", filename)
                instance = IPInstance(filename)
                ipgbs_result = IPGBs.groebner_basis(instance)
                fourti2_result = GBTools.tovector(IPGBs.FourTi2.groebner(instance))
                test_buchberger_correctness(ipgbs_result, fourti2_result, instance)
            end
        end
    end

end
