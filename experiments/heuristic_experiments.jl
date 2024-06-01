using IPGBs
using MIPMatrixTools.IPInstances
using MIPMatrixTools.MatrixTools

using Graphs
using JuMP
using Random

function dominating_set_instance(n :: Int, m :: Int)
    G = SimpleGraph(n)
    edges = 0
    while edges < m
        u = rand(1:n)
        v = rand(1:n)
        if u != v && !has_edge(G, u, v)
            add_edge!(G, u, v)
            edges += 1
        end
    end
    model = Model()
    @variable(model, x[1:n], Bin)
    @constraint(model, [i in 1:n], x[i] + sum(x[j] for j in neighbors(G, i)) >= 1)
    @objective(model, Min, sum(x))
    ip = IPInstance(model)
    return ip
end

function experiment(;p = 0.1, reps = 30, time_limit = 10.0, gb_size_limit = 100000)
    Random.seed!(0)
    println("n m IPGBs OPT")
    for n in 10:10:100
        for _ in 1:reps
            m = round(Int, p * binomial(n, 2))
            ip = dominating_set_instance(n, m)
            init_sol = ones(Int, n)
            full_sol = MatrixTools.lift_partial_solution(init_sol, ip.b, ip.A)
            gb_heuristic!(
                full_sol, ip, time_limit = time_limit, gb_size_limit = gb_size_limit
            )
            val = round(Int, ip.C[1, :]' * full_sol)
            _, opt_val, _ = IPInstances.solve(ip)
            println(n, " ", m, " ", val, " ", opt_val)
        end
    end
end
