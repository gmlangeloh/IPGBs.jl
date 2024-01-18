"""
Incomplete implementation of a Zero-suppressed decision diagram (ZDD).

In theory, ZDDs could be an efficient way to implement the reducer search, but implementing
this data structure well is a bit tricky. I intend to complete this implementation later.
"""
module ZDDs

export ZDD

using DataStructures
using Graphs
using IPGBs.FastBitSets

struct ZDDNode
    index :: Int
    low :: Union{ZDDNode, Nothing}
    high :: Union{ZDDNode, Nothing}
end

const TOP = ZDDNode(0, nothing, nothing)
const BOTTOM = ZDDNode(-1, nothing, nothing)

is_top(node :: ZDDNode) = node.index == 0
is_bottom(node :: ZDDNode) = node.index == -1

struct ZDD
    nodes :: Vector{ZDDNode}
    root :: ZDDNode
end

const TOP_ZDD = ZDD([TOP], TOP)
const BOTTOM_ZDD = ZDD([BOTTOM], BOTTOM)

function to_digraph(zdd :: ZDD) :: DiGraph
    #Keep track of the node labels / indices for printing later
    node_labels = Int[]
    visited_nodes = Set{ZDDNode}()
    D = DiGraph()
    #Traverse the ZDD in a depth-first manner, adding visited nodes and edges to D
    s = Stack{ZDDNode}()
    push!(s, zdd.root)

    while !isempty(s)
        node = pop!(s)
        if !visited_nodes[node]
            visited_nodes[node] = true
            push!(node_labels, node.index)
            if !is_top(node)
                push!(s, node.low)
                push!(s, node.high)
                add_edge!(D, node.index, node.low.index)
                add_edge!(D, node.index, node.high.index)
            end
        end
    end
    return D
end

function _partition_by_index(sets :: Vector{FastBitSet}, index :: Int)
    low = Vector{FastBitSet}()
    high = Vector{FastBitSet}()
    for set in sets
        if set[index]
            push!(high, set)
        else
            push!(low, set)
        end
    end
    return low, high
end

#TODO: Still need to make sure the ZDD is reduced. There are currently redundant nodes
function _build_ZDD(
    sets :: Vector{FastBitSet},
    n :: Int,
    i :: Int = 1
) :: ZDD
    if isempty(sets)
        return BOTTOM_ZDD
    elseif i > n
        return TOP_ZDD
    end
    #Find smallest j >= i such that some set in sets contains j
    j = i
    found = false
    while !found && j <= n
        for set in sets
            if set[j]
                found = true
                break
            end
        end
        if !found
            j += 1
        end
    end
    if j > n
        return TOP_ZDD
    end
    #Take the sets that don't contain j, build ZDD recursively
    #Take the sets that do contain j, build ZDD recursively
    low, high = _partition_by_index(sets, j)
    low_zdd = _build_ZDD(low, n, j + 1)
    high_zdd = _build_ZDD(high, n, j + 1)
    #Connect node j to the two recursively built ZDDs
    j_node = ZDDNode(j, low_zdd.root, high_zdd.root)
    all_nodes = vcat(low_zdd.nodes, high_zdd.nodes, [j_node])
    #Remove redundant nodes
    return ZDD(all_nodes, j_node)
end

function ZDD(sets :: Vector{FastBitSet})
    n = 0
    if !isempty(sets)
        n = length(sets[1])
    end
    return _build_ZDD(sets, n)
end

end
