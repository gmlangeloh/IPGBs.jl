using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.FourTi2
using IPGBs.GBAlgorithms
using IPGBs.GBTools
using IPGBs.MultiObjectiveAlgorithms
using IPGBs.MultiObjectiveTools
using MultiObjectiveInstances
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

    include("./walkback_tests.jl")

    @testset "Guessing initial solutions" begin
        for filename in readdir("test_instances", join=true)
            if endswith(filename, ".mps")
                println("Initial solution test for ", filename)
                instance = IPInstance(filename)
                solution = IPInstances.guess_initial_solution(instance)
                println("Guessed feasible solution? ",
                    @test is_feasible_solution(instance, solution))
                println()
            end
        end
    end

    @testset "Optimize combinatorial optimization instances" begin
        for filename in readdir("test_instances", join=true)
            if endswith(filename, ".mps")
                println("Optimization test for ", filename)
                instance = IPInstance(filename)
                #Optimize with IPGBs and with a traditional IP solver, then compare
                _, ipgbs_value, ipgbs_status = optimize(instance)
                solver_solution, solver_value, solver_status = IPInstances.solve(instance)
                println("Optimal value with no initial solution? ",
                    @test compare_to_solver(ipgbs_value, ipgbs_status, solver_value, solver_status))
                init_solution = IPInstances.guess_initial_solution(instance)
                _, init_value, init_status = optimize(instance, solution=init_solution)
                println("Optimal value starting from arbitrary solution? ",
                    @test compare_to_solver(init_value, init_status, solver_value, solver_status))
                _, opt_value, opt_status = optimize(instance, solution=solver_solution)
                println("Optimal value starting from the optimal solution? ",
                    @test compare_to_solver(opt_value, opt_status, solver_value, solver_status))
                println()
            end
        end
    end

    @testset "MOIP tests" begin
        for filename in readdir("test_instances_moip", join=true)
            if endswith(filename, ".dat")
                instance = MultiObjectiveInstances.read_from_file(filename)
                initial_solution = MultiObjectiveInstances.Knapsack.knapsack_initial(instance)
                for solver in ["4ti2", "IPGBs"]
                    println("MOIP test for ", filename, " with solver ", solver)
                    Xs, Ys, _ = moip_gb_solve(instance, initial_solution, solver=solver)
                    println("Is the efficient set feasible? ",
                        @test is_efficient_set_feasible(Xs, instance.A, instance.b)
                    )
                    println("Is the Pareto set nondominated? ",
                        @test is_nondominated_set(Ys, maximization=false)
                    )
                    if size(instance.A, 2) < 15
                        #Compare the results to the very slow enumeration algorithm
                        X2s, Y2s = moip_walkback(instance)
                        println("Same nondominated set as enumeration? ",
                            @test Ys == Y2s
                        )
                        println("Same efficient set as enumeration? ",
                            @test Set(Xs) == Set(X2s)
                        )
                    end
                    println()
                end
            end
        end
    end
end
