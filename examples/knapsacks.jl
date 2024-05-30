# In this example, we explore GBs of unbounded / binary knapsacks
# We also generate some figures of their digraphs to illustrate the results
# The reason this is not a notebook is that plotting interacts weirdly with VSCode
# notebooks in Julia. It appears to be a bug in Jupyter or something
# https://github.com/JuliaPlots/Plots.jl/issues/4827

using IPGBs
using IPGBs.Markov
using MIPMatrixTools.IPInstances

using Graphs
using MetaGraphs
using GraphMakie
using Makie
using CairoMakie
using NetworkLayout
using LayeredLayouts

function basic_knapsack(binary = false)
    A = [15 18 19 17 20]
    b = [45]
    C = [1 1 1 1 1]
    u = [binary ? 1 : nothing for i in 1:5]
    return IPInstance(A, b, C, u)
end

function feasible_set(A, b, binary = false)
    #TODO: Ideally I should use the Walkback method
    #Currently I assume the instance is binary_knapsack
    feas = Vector{Int}[[0, 0, 0, 0, 0, 45]]
    #Compute all solutions with a single non-zero variable
    for i in 1:5
        for val in 1:3
            if binary && val > 1
                break
            end
            if val * A[1, i] <= b[1]
                sol = zeros(Int, 6)
                sol[i] = val
                sol[6] = b[1] - val * A[1, i]
                push!(feas, sol)
            end
        end
    end
    #Compute all solutions with two non-zero variables
    for i in 1:5
        for j in 1:5
            if i > j
                if A[1, i] + A[1, j] <= b[1]
                    sol = zeros(Int, 6)
                    sol[i] = 1
                    sol[j] = 1
                    sol[6] = b[1] - A[1, i] - A[1, j]
                    push!(feas, sol)
                end
            end
        end
    end
    #There can be no other feasible solutions for these instances
    return feas
end

function gb_skeleton(solutions, gb)
    skeleton = SimpleDiGraph(length(solutions))
    #Add vertex labels corresponding to the solutions
    #for i in 1:length(solutions)
    #    set_prop!(skeleton, i, :label, solutions[i])
    #end
    #Add edges corresponding to the elements in the GB
    edge_indices = Int[]
    for i in eachindex(solutions)
        for j in eachindex(solutions)
            if i == j
                continue
            end
            k = findfirst(isequal(solutions[i] - solutions[j]), gb)
            if !isnothing(k)
                add_edge!(skeleton, i, j)
                push!(edge_indices, k)
            end
        end
    end
    return skeleton, edge_indices
end

function lift_binary_solution(solution)
    original_vars = solution[1:5]
    return [solution; -original_vars]
end

function diverse_color(i, n)
    r = sin(2pi*i/n + 0.0)
    g = sin(2pi*i/n + 2pi/3)
    b = sin(2pi*i/n + 4pi/3)
    return Makie.RGBA(r, g, b)
end

function optimal_solution(skeleton)
    return findall(v -> length(outneighbors(skeleton, v)) == 0, vertices(skeleton))[1]
end

function is_binary(x)
    return all(x[i] == 1 || x[i] == 0 for i in 1:5)
end

function plot_int_skeleton()
    int_knapsack = basic_knapsack(false)
    int_solutions = feasible_set(int_knapsack.A, int_knapsack.b, false)
    int_gb = groebner_basis(int_knapsack)
    int_skeleton, edge_indices = gb_skeleton(int_solutions, int_gb)
    opt_index = optimal_solution(int_skeleton)
    println(opt_index, ": ", int_solutions[opt_index])
    edge_colors = [diverse_color(i, length(int_gb)) for i in edge_indices]
    graphplot(
        int_skeleton,
        nlabels=[string(i) for i in 1:length(int_solutions)],
        node_color=[is_binary(v) for v in int_solutions],
        edge_color=edge_colors,
        layout=Shell()
    )
end

function plot_bin_skeleton()
    bin_knapsack = basic_knapsack(true)
    bin_solutions = feasible_set(bin_knapsack.A, bin_knapsack.b, true)
    bin_gb = groebner_basis(bin_knapsack)
    bin_skeleton, edge_indices = gb_skeleton(lift_binary_solution.(bin_solutions), bin_gb)
    opt_index = optimal_solution(bin_skeleton)
    println(opt_index, ": ", bin_solutions[opt_index])
    edge_colors = [diverse_color(i, length(bin_gb)) for i in edge_indices]
    graphplot(
        bin_skeleton,
        nlabels=[string(i) for i in 1:length(bin_solutions)],
        edge_color=edge_colors,
        layout=Shell()
    )
end
