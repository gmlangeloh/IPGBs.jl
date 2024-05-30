module FeasibleGraphs

using IPGBs
using IPGBs.Walkback
using MIPMatrixTools.IPInstances

using CairoMakie
using Graphs
using GraphMakie
using LayeredLayouts
using Makie
using MetaGraphs
using Random

export feasible_graph, plot_feasible_graph

function feasible_graph(instance :: IPInstance)
    gb = groebner_basis(instance)
    return feasible_graph(instance, gb)
end

function feasible_graph(instance :: IPInstance, gb :: Vector{Vector{Int}})
    feasible_solutions = enumerate_solutions(instance)
    graph = MetaDiGraph(length(feasible_solutions))
    edge_indices = Int[]
    vector_sols = [ sol for sol in feasible_solutions ]
    for i in eachindex(vector_sols)
        set_prop!(graph, i, :label, vector_sols[i])
    end
    for i in eachindex(vector_sols)
        for j in eachindex(vector_sols)
            if i == j
                continue
            end
            k = findfirst(isequal(vector_sols[i] - vector_sols[j]), gb)
            if !isnothing(k)
                add_edge!(graph, i, j)
                push!(edge_indices, k)
            end
        end
    end
    return graph, edge_indices
end

function plot_feasible_graph(instance :: IPInstance)
    graph, edge_indices = feasible_graph(instance)
    plot_feasible_graph(graph, edge_indices)
end

function plot_feasible_graph(
    instance :: IPInstance,
    arc_set :: Vector{Vector{Int}}
)
    graph, edge_indices = feasible_graph(instance, arc_set)
    plot_feasible_graph(graph, edge_indices)
end

function random_colors(n)
    Random.seed!(1)
    colors = []
    for i in 1:n
        r = rand(0:255)
        g = rand(0:255)
        b = rand(0:255)
        c = Makie.RGBA{Makie.N0f8}(r / 255, g / 255, b / 255, 255 / 255)
        push!(colors, c)
    end
    return colors
end

function dag_layout(g :: AbstractGraph)
   xs, ys, _ = solve_positions(Zarate(), g)
   return Point.(zip(xs, ys))
end

function plot_feasible_graph(
    graph :: AbstractGraph,
    edge_indices :: Vector{Int};
    full_labels :: Bool = false
)
    colors = random_colors(length(edge_indices))
    labels = [string(i) for i in 1:nv(graph)]
    if full_labels
        labels = [string(get_prop(graph, i, :label)) for i in 1:nv(graph)]
    end
    f, ax, p = graphplot(
        graph,
        nlabels=labels,
        edge_color=[colors[i] for i in edge_indices],
        layout=dag_layout,
    )
    hidedecorations!(ax); hidespines!(ax)
    f
end

end
