using Test

include("./test_functions.jl")

@testset "IPGBs.jl" begin
    #Tests for Buchberger implementation
    for s in [ false, true ]
        for n in [5, 10, 15, 20, 25]
            s_name = s ? "GradedBinomial" : "Binomial"
            println("Buchberger test for ", s, " structure, n = ", n)
            gb, fourti2gb = test_buchberger(n, implicit_representation=s)
            @test length(gb) == length(fourti2gb)
            @test IPGBs.GBTools.isequal(gb, fourti2gb)
            println()
        end
    end

    for n in [5, 10, 15, 20]
        println("Signature test for Binomial structure, n = ", n)
        gb, fourti2gb, bgb = test_siggb(n)
        @test IPGBs.GBTools.isincluded(bgb, gb)
        println()
    end

    #Tests for the signature-based algorithms: compare to Buchberger
    #for s in [ Binomial, GradedBinomial ]
    #    for n in [5, 10, 15, 20, 25]
    #        println("Signature test for ", s, " structure, n = ", n)
    #        gb, fourti2gb, bgb = test_siggb(n, structure=s)
    #        @test length(gb) == length(bgb)
    #        @test IPGBs.GBTools.isequal(gb, bgb)
    #        println()
    #    end
    #end
end
