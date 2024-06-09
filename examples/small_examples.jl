using IPGBs
using IPGBs.Markov
using IPGBs.FeasibleGraphs

using IPGBs.IPInstances

using Graphs
using JuMP
using Random

Random.seed!(0)

function find_small_example(
    n :: Int, m :: Int, max_solutions :: Int, min_sinks :: Int, binary=true
)
    ip = nothing
    mb = nothing
    g = nothing
    example_found = false
    while !example_found
        #Generate random binary knapsack with correlation between values and weights
        C = rand(2:10, 1, n)
        A = rand(2:10, m, n)
        b = [round(Int, sum(A[i, :]) / 2) for i in 1:m]
        for j in 1:n
            A[j] = rand(C[1, j] - 5:C[1, j] + 5)
        end
        u = [binary ? 1 : nothing for _ in 1:n]
        ip = IPInstance(A, b, C, u)
        #Compute its feasibility graph with respect to the Markov basis
        mb = markov_basis(ip)
        g, _ = feasible_graph(ip, mb)
        #Check whether the knapsack has at most max_solutions feasible points
        #and the minimum number of sinks. Keep generating new knapsacks until this condition
        #is satisfied.
        if nv(g) <= max_solutions && count(d == 0 for d in outdegree(g)) >= min_sinks
            example_found = true
        end
    end
    println("Found graph with:")
    println("Markov basis size: ", length(mb))
    println("Number of vertices: ", nv(g))
    println("Number of sinks: ", count(d == 0 for d in outdegree(g)))
    return ip
end

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

function evaluate_dominating_set(n :: Int, m :: Int)
    ip = dominating_set_instance(n, m)
    mb = markov_basis(ip)
    g, _ = feasible_graph(ip, mb)
    return ip
end
