using Test

include("./test_functions.jl")

@testset "IPGBs.jl" begin
    #Tests for Buchberger implementation
    for s in [ false, true ]
        for n in [5, 10, 15, 20, 25]
            s_name = s ? "GradedBinomial" : "Binomial"
            println("Buchberger test for ", s_name, " structure, n = ", n)
            gb, fourti2gb = test_buchberger(n, implicit_representation=s)
            @test length(gb) == length(fourti2gb)
            @test IPGBs.GBTools.isequal(gb, fourti2gb)
            println()
        end
    end

    for order in [:ltpot, :pot, :top]
        for n in [5, 10, 15, 20]
            println("Signature test for Binomial structure, n = ", n)
            gb, fourti2gb, bgb = test_siggb(n, module_order=order)
            @test IPGBs.GBTools.isincluded(bgb, gb)
            println()
        end
    end
end
