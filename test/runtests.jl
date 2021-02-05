using Test

include("./test_functions.jl")

@testset "Buchberger: is_truncated_groebner_basis" begin
    for s in [ false, true ]
        for n in [5, 10, 15, 20, 25, 30]
            s_name = s ? "GradedBinomial" : "Binomial"
            println("Buchberger test for ", s_name, " structure, n = ", n)
            gb, A, b, C, u = run_algorithm(n, use_signatures=false, implicit_representation=s,
                                           quiet=true)
            bs = BinomialSet(gb, C)
            @test is_truncated_groebner_basis(bs, A, b, u)
            println()
        end
    end
end

@testset "Bucheberger: is_groebner_basis" begin
    for s in [ false, true ]
        for n in [5, 10, 15, 20]
            s_name = s ? "GradedBinomial" : "Binomial"
            println("Buchberger test for ", s_name, " structure, n = ", n)
            gb, A, b, C, u = run_algorithm(n, use_signatures=false, implicit_representation=s,
                                           quiet=true, truncate=false)
            bs = BinomialSet(gb, C)
            @test is_groebner_basis(bs)
            println()
        end
    end
end

@testset "Signature: is_truncated_groebner_basis" begin
    for order in [:ltpot, :pot, :top]
        for n in [5, 10, 15, 20]
            println("Signature test for Binomial structure, n = ", n)
            gb, A, b, C, u = run_algorithm(n, module_order=order, quiet=true)
            bs = BinomialSet(gb, C)
            @test is_truncated_groebner_basis(bs, A, b, u)
            println()
        end
    end
end

@testset "Signature: is_groebner_basis" begin
    for order in [:ltpot, :pot, :top]
        for n in [5, 10, 15, 20]
            println("Signature test for Binomial structure, n = ", n)
            gb, A, b, C, u = run_algorithm(n, module_order=order, quiet=true, truncate=false)
            bs = BinomialSet(gb, C)
            @test is_groebner_basis(bs)
            println()
        end
    end
end

@testset "Buchberger: compare with 4ti2" begin
    println("Starting Buchberger tests")
    println("------------------------------------------------")
    for s in [ false, true ]
        for n in [5, 10, 15, 20, 25, 30]
            s_name = s ? "GradedBinomial" : "Binomial"
            println("Buchberger test for ", s_name, " structure, n = ", n)
            gb, fourti2gb = test_buchberger(n, implicit_representation=s, quiet=true)
            @test length(gb) == length(fourti2gb)
            @test IPGBs.GBTools.isequal(gb, fourti2gb)
            println()
        end
    end
end

@testset "Signature: compare with Buchberger" begin
    for order in [:ltpot, :pot, :top]
        for n in [5, 10, 15, 20]
            println("Signature test for Binomial structure, n = ", n)
            gb, fourti2gb, bgb = test_siggb(n, module_order=order, quiet=true)
            @test IPGBs.GBTools.isincluded(bgb, gb)
            println()
        end
    end
end
