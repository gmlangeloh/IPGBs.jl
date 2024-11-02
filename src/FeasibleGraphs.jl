module FeasibleGraphs

using IPGBs
using IPGBs.Orders
using IPGBs.Walkback
using IPGBs.IPInstances

using GLMakie
using Graphs
using GraphMakie
using LayeredLayouts
using Makie
using Makie.Colors
using MetaGraphs
using Random

export feasible_graph, plot_feasible_graph

function feasible_graph(instance :: IPInstance; max_vertices :: Int = 0)
    gb = groebner_basis(instance)
    return feasible_graph(instance, gb, max_vertices = max_vertices)
end

function feasible_graph(
    instance :: IPInstance,
    arc_set :: Vector{Vector{Int}};
    max_vertices :: Int = 0
)
    feasible_solutions = enumerate_solutions(instance, max_solutions=max_vertices)
    #Orient arc_set according to instance
    oriented_arcs = Vector{Vector{Int}}()
    order = MonomialOrder(instance)
    for arc in arc_set
        if Orders.is_inverted_generic(order, arc)
            push!(oriented_arcs, -arc)
        else
            push!(oriented_arcs, arc)
        end
    end
    return feasible_graph(feasible_solutions, oriented_arcs)
end

function feasible_graph(
    feasible_solutions :: Set{Vector{Int}},
    arc_set :: Vector{Vector{Int}}
)
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
            k = findfirst(isequal(vector_sols[i] - vector_sols[j]), arc_set)
            if !isnothing(k)
                add_edge!(graph, i, j)
                push!(edge_indices, k)
            end
        end
    end
    return graph, edge_indices
end

false_predicate(x) = false

function plot_feasible_graph(
    instance :: IPInstance;
    max_vertices :: Int = 0,
    full_labels :: Bool = false,
    predicate = false_predicate
)
    graph, edge_indices = feasible_graph(instance, max_vertices = max_vertices)
    plot_feasible_graph(graph, edge_indices, full_labels = full_labels, predicate=predicate)
end

function plot_feasible_graph(
    instance :: IPInstance,
    arc_set :: Vector{Vector{Int}};
    max_vertices :: Int = 0,
    full_labels :: Bool = false,
    predicate = false_predicate
)
    graph, edge_indices = feasible_graph(instance, arc_set, max_vertices = max_vertices)
    plot_feasible_graph(graph, edge_indices, full_labels = full_labels, predicate=predicate)
end

function plot_feasible_graph(
    feasible_solutions :: Vector{Vector{Int}},
    arc_set :: Vector{Vector{Int}};
    full_labels :: Bool = false,
    predicate = false_predicate
)
    graph, edge_indices = feasible_graph(feasible_solutions, arc_set)
    plot_feasible_graph(graph, edge_indices, full_labels = full_labels, predicate=predicate)
end

function dag_layout(g :: AbstractGraph)
   xs, ys, _ = solve_positions(Zarate(), g)
   return Point.(zip(xs, ys))
end

function plot_feasible_graph(
    graph :: AbstractGraph,
    edge_indices :: Vector{Int};
    full_labels :: Bool = false,
    predicate = false_predicate
)
    #colors = random_colors(maximum(edge_indices))
    hsv_colors = HSV.(range(0, 360, length(edge_indices)), 10, 10)
    labels = [string(i) for i in 1:nv(graph)]
    if full_labels
        labels = [string(get_prop(graph, i, :label)) for i in 1:nv(graph)]
    end
    f, ax, _ = graphplot(
        graph,
        nlabels=labels,
        node_color=[predicate(get_prop(graph, i, :label)) ? :red : :black for i in 1:nv(graph)],
        arrow_size=25,
        edge_color=[hsv_colors[i] for i in edge_indices],
        layout=dag_layout,
    )
    hidedecorations!(ax); hidespines!(ax)
    f
end

end
