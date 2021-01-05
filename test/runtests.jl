using Test

include("./test_functions.jl")

@testset "IPGBs.jl" begin
    for s in [ Binomial, GradedBinomial ]
        for n in [5, 10, 15, 20, 25]
            println("Buchberger test for ", s, " structure, n = ", n)
            gb, fourti2gb = test_buchberger(n, structure=s)
            @test length(gb) == length(fourti2gb)
            @test IPGBs.GBTools.isequal(gb, fourti2gb)
            println()
        end
    end

    #TODO Add tests for siggb
end
