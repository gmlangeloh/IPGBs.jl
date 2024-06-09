using IPGBs
using IPGBs.IPInstances
using IPGBs.MatrixTools

using Graphs
using JuMP
using Random

function dominating_set_instance(n :: Int, m :: Int, filename :: String = "")
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
    if filename != ""
        write_to_file(model, filename)
    end
    ip = IPInstance(model)
    return ip
end

function experiment(;p = 0.1, reps = 5, time_limit = 10.0, gb_size_limit = 100000)
    Random.seed!(0)
    println("n m IPGBs OPT IPGBs_time Solver_time")
    for n in 10:10:100
        for r in 1:reps
            m = round(Int, p * binomial(n, 2))
            filename = "dominating_set/dominating_set_$(n)_$(m)_$(r).mps"
            ip = dominating_set_instance(n, m, filename)
            init_sol = ones(Int, n)
            full_sol = MatrixTools.lift_partial_solution(init_sol, ip.b, ip.A)
            ipgbs_stats = @timed gb_heuristic!(
                full_sol, ip, time_limit = time_limit, gb_size_limit = gb_size_limit
            )
            val = round(Int, ip.C[1, :]' * full_sol)
            solver_stats = @timed IPInstances.solve(ip)
            opt_val = solver_stats.value[2]
            println(
                n, " ", m, " ", val, " ", opt_val, " ",
                ipgbs_stats.time, " ", solver_stats.time
            )
        end
    end
end

function larger_experiment(;p = 0.1, reps = 5, time_limit = 10.0, gb_size_limit = 100000)
    Random.seed!(0)
    println("n m IPGBs OPT IPGBs_time Solver_time")
    for n in 100:100:1000
        for r in 1:reps
            m = round(Int, p * binomial(n, 2))
            filename = "dominating_set/dominating_set_$(n)_$(m)_$(r).mps"
            ip = dominating_set_instance(n, m, filename)
            init_sol = ones(Int, n)
            full_sol = MatrixTools.lift_partial_solution(init_sol, ip.b, ip.A)
            ipgbs_stats = @timed gb_heuristic!(
                full_sol, ip, time_limit = time_limit, gb_size_limit = gb_size_limit
            )
            val = round(Int, ip.C[1, :]' * full_sol)
            solver_stats = @timed IPInstances.solve(ip)
            opt_val = solver_stats.value[2]
            println(
                n, " ", m, " ", val, " ", opt_val, " ",
                ipgbs_stats.time, " ", solver_stats.time
            )
        end
    end
end

function dominated_set_markov(n :: Int, ip :: IPInstance)
    markov = Vector{Int}[]
    for i in 1:n
        m = zeros(Int, n)
        m[i] = 1
        rhs = zeros(Int, ip.m)
        lifted_m = MatrixTools.lift_partial_solution(m, rhs, ip.A)
        push!(markov, lifted_m)
    end
    return markov
end

function dominating_set_experiment(filename; time_limit = 10.0, gb_size_limit = 100000)
    ip = IPInstance(filename)
    n = round(Int, ip.m / 2)
    init_sol = ones(Int, n)
    full_sol = MatrixTools.lift_partial_solution(init_sol, ip.b, ip.A)
    mb = dominated_set_markov(n, ip)
    ipgbs_stats = @timed gb_heuristic!(
        full_sol, ip, mb, time_limit = time_limit, gb_size_limit = gb_size_limit
    )
    val = round(Int, ip.C[1, :]' * full_sol)
    println(ipgbs_stats.time, " ", val)
end
