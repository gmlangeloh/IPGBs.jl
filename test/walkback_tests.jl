using IPGBs.Walkback
using Test

function enumeration_tests()
    @testset "Walkback enumeration" begin
        A = [15 18 19 17 20]
        b = [45]
        C = [1 1 1 1 1]
        u = [nothing for _ in 1:5]
        #Test unbounded knapsack
        ip = IPInstance(A, b, C, u)
        sols = enumerate_solutions(ip)
        println("Correct feasible solutions in integer instance? ",
            @test length(sols) == 22
        )

        #Test binary knapsack
        u_bin = [1 for _ in 1:5]
        bin_ip = IPInstance(A, b, C, u_bin)
        bin_sols = enumerate_solutions(bin_ip)
        println("Correct feasible solutions in binary instance? ",
            @test length(bin_sols) == 16
        )

        println()
    end
end

enumeration_tests()
