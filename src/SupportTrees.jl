"""
Reimplementation of 4ti2's FilterReduction tree structure used to store
binomial sets for efficient reductions. Called SupportTrees in reference
to Roune and Stillman (2012).
"""
module SupportTrees
export SupportTree, find_reducer, support_tree, addbinomial!,
    removebinomial!, enumerate_reducers, find_reducer_iter

using IPGBs.GBElements
using IPGBs.Statistics

mutable struct TreeStats <: GBStats
    reduction_steps :: Int #Number of calls to find_reducer (recursive)
    reducers_checked :: Int #Number of candidate reducers / calls to GBElements.reduces

    TreeStats() = new(0, 0)
end

"""
Each binomial in the GB is inserted in some node of the tree.

A filter is a list of indices of the binomials such that their coordinates
are positive - that is, the list of indices of variables appearing in its
leading term.

A binomial is inserted in a node of the tree with the same filter as itself.
All binomials with the same filter are inserted in the same node. If the
tree corresponds to a GB (or a partial GB in some way) each node should
contain at most one binomial in its binomial_list.

Each child node is labeled by its (next) filter index.
"""
struct SupportNode{T <: AbstractVector{Int}}
    children :: Vector{Tuple{Int, SupportNode{T}}}
    binomial_list :: Vector{T}
    filter :: Vector{Int}
    SupportNode{T}() where {T <: AbstractVector} = new(Tuple{Int, SupportNode{T}}[], T[], Int[])
end

function Base.show(
    io :: IO,
    node :: SupportNode{T}
) where T
    num_binomials = length(node.binomial_list)
    num_children = length(node.children)
    println(io, "Node filter: ", node.filter)
    println(io, "Node with ", num_binomials, " elements, ",
            num_children, " children")
    for binomial in node.binomial_list
        println(io, binomial)
    end
    println(io, "End node.")
    for (index, child) in node.children
        println(io)
        println(io, "Child index: ", index)
        print(io, child)
    end
end

"""
The structure from 4ti2's FilterReduction for efficient reduction of binomials
in toric ideals.

At each level `i` of the tree, there exists a node labeled `j` if there is some
binomial in the tree whose filter's i-th element is j. Each binomial is stored
at level `n`, with `n` the size of its filter.

If `fullfilter` is true, then both positive and negative entries of the
binomials are considered to be in the filter for reduction. This is useful to
implement the specialized truncated GB algorithm described in Thomas and
Weismantel, Section 3.
"""
mutable struct SupportTree{T <: AbstractVector{Int}}
    root :: SupportNode{T}
    size :: Int
    depth :: Int #This is useful to know how balanced the tree is
    fullfilter :: Bool
    stats :: TreeStats

    #For efficiency of the find_reducer_iter function
    node_stack :: Vector{Tuple{Int, SupportNode{T}}}

    SupportTree{T}(fullfilter) where {T <: AbstractVector{Int}} = begin
        new(SupportNode{T}(), 0, 1, fullfilter, TreeStats(), Vector(undef, 2))
    end
end

function Base.show(
    io :: IO,
    tree :: SupportTree{T}
) where T
    println(io, "Support tree with ", tree.size, " elements")
    print(io, tree.root)
end

"""
Puts references to the elements of `gb` in a reduction tree, so that one may
easily find whether a given binomial can be reduced by `gb`.
"""
function support_tree(
    gb :: S;
    fullfilter :: Bool = false
) :: SupportTree{T} where { T <: AbstractVector{Int}, S <: AbstractVector{T}}
    tree = SupportTree{T}(fullfilter)
    for i in eachindex(gb)
        addbinomial!(tree, gb[i])
    end
    return tree
end

"""
Adds a (reference to a) binomial to a SupportTree.

The idea is adding the binomial to a node corresponding precisely to its
filter. At each iteration, take the i-th element of the filter of `binomial`.
If there is already a node in the tree at level i labeled by the i-th index
in the filter of binomial, move to the subtree defined by that node.
Otherwise, create such a node.
"""
function addbinomial!(
    tree :: SupportTree{T},
    binomial :: T
) where {T <: AbstractVector{Int}}
    current = tree.root
    depth = 1
    binomial_filter = GBElements.filter(binomial, fullfilter=tree.fullfilter)
    #Search for a path in the tree labeled by the indices in the filter
    #Create any missing nodes in this path
    for i in binomial_filter
        j = 1
        while j <= length(current.children) && current.children[j][1] != i
            j += 1
        end
        if j <= length(current.children) #Next node in the path already exists
            current = current.children[j][2]
        else #Create new node in the tree
            newnode = SupportNode{T}()
            push!(current.children, (i, newnode))
            current = last(current.children)[2]
            #break #I think this is enough, it leaves the loop anyway
        end
        depth += 1
    end
    push!(current.binomial_list, binomial)
    if isempty(current.filter) #Set current node's filter if it wasn't already
        for i in binomial_filter #current.filter = f
            push!(current.filter, i)
        end
    end
    tree.depth = max(tree.depth, depth)
    resize!(tree.node_stack, tree.depth + 1)
    tree.size += 1
end

"""
Removes `binomial` from `tree` by doing a lookup operation and then just
removes it from the corresponding binomial_list.

If `binomial` is not in `tree`, nothing will happen.
"""
function removebinomial!(
    tree :: SupportTree{T},
    binomial :: T
) where {T <: AbstractVector{Int}}
    current = tree.root
    binomial_filter = GBElements.filter(binomial, fullfilter=tree.fullfilter)
    for i in binomial_filter
        j = 1
        while j <= length(current.children) && current.children[j][1] != i
            j += 1
        end
        if j <= length(current.children)
            current = current.children[j][2]
        end
    end
    for i in 1:length(current.binomial_list)
        if current.binomial_list[i] === binomial
            deleteat!(current.binomial_list, i)
            return
        end
    end
end

"""
Returns a vector containing all reducers of `g` in `gb`, using `tree` to find them efficiently.
"""
function enumerate_reducers(
    g :: P,
    gb :: S,
    tree :: SupportTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Vector{T} where {P <: AbstractVector{Int}, T <: AbstractVector{Int}, S <: AbstractVector{T}}
    reducers = T[]
    is_singular = Ref{Bool}(false) #Useless for now
    enumerate_reducers!(
        reducers, g, gb, tree.root, fullfilter=tree.fullfilter,
        skipbinomial=skipbinomial, negative=negative,
        is_singular=is_singular
    )
    return reducers
end

"""
Walks in the tree recursively pushing any reducers found to `reducers`.
"""
function enumerate_reducers!(
    reducers :: Vector{T},
    g :: P,
    gb :: S,
    node :: SupportNode{T};
    fullfilter :: Bool = false,
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false,
    is_singular :: Ref{Bool} = Ref(false)
) where {P <: AbstractVector{Int}, T <: AbstractVector{Int}, S <: AbstractVector{T}}
    for (i, child) in node.children
        if g[i] > 0 || (negative && g[i] < 0) || (fullfilter && g[i] != 0)
            #Look for reducers recursively in the children of this node
            enumerate_reducers!(
                reducers, g, gb, child, fullfilter=fullfilter,
                skipbinomial=skipbinomial, negative=negative,
                is_singular=is_singular
            )
        end
    end
    #Check if some binomial in this node works as reducer
    for reducer in node.binomial_list
        # If the reducer is the same element of the GB that was passed as
        #parameter to be skipped, skip it.
        # This is useful in inter-reductions, where the element should not be
        #used to reduce itself.
        if !isnothing(skipbinomial) && (reducer === skipbinomial || isequal(reducer, skipbinomial))
            continue
        end
        if GBElements.reduces(
            g, node.filter, reducer, gb, fullfilter=fullfilter, negative=negative,
            is_singular=is_singular
        )
            push!(reducers, reducer)
        end
    end
end

"""
Finds a reducer of `g` in `tree` by searching the tree iteratively, instead of
recursively. It serves as an alternative to find_reducer, although the performance
is currently nearly equivalent.
"""
function find_reducer_iter(
    g :: T,
    gb :: S,
    tree :: SupportTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false,
    is_singular :: Ref{Bool} = Ref(false)
) :: Union{T, Nothing} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    node_stack = tree.node_stack
    stack_index = 1
    node_stack[stack_index] = (1, tree.root)
    #push!(node_stack, (1, tree.root))
    while stack_index >= 1 #!isempty(node_stack)
        child_index, current_node = node_stack[stack_index] #pop!(node_stack)
        stack_index -= 1
        if child_index <= length(current_node.children)
            stack_index += 1
            node_stack[stack_index] = (child_index + 1, current_node)
            #push!(node_stack, (child_index + 1, current_node))
            i, child = current_node.children[child_index]
            if g[i] > 0 || negative && g[i] < 0 || (tree.fullfilter && g[i] != 0)
                #push!(node_stack, (1, child))
                stack_index += 1
                node_stack[stack_index] = (1, child)
            end
        else #Check whether this node contains a reducer
            for reducer in current_node.binomial_list
                if !isnothing(skipbinomial) && reducer === skipbinomial
                    continue
                end
                if GBElements.reduces(
                    g, current_node.filter, reducer, gb, fullfilter=tree.fullfilter,
                    negative=negative, is_singular=is_singular
                )
                    return reducer
                end
            end
        end
    end
    return nothing
end

"""
Checks whether `g` is reducible by some element of `gb` and, if so, returns
a reducer. Otherwise, returns `nothing`.
"""
function find_reducer(
    g :: T,
    gb :: S,
    tree :: SupportTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false,
    is_singular :: Ref{Bool} = Ref(false)
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    return find_reducer(
        g, gb, tree.root, tree.stats, fullfilter=tree.fullfilter,
        skipbinomial=skipbinomial, negative=negative, is_singular=is_singular
    )
end

"""
Checks whether `g` is reducible by some element of `gb` and, if so, returns
a reducer. Otherwise, returns `nothing`.

This is done by looking up a node in this subtree with a subfilter of g's
filter and then checking directly whether some binomial in the node reduces
g.

The lookup operation traverses the tree looking for nodes whose labels are
contained in the filter of `g`.
"""
function find_reducer(
    g :: T,
    gb :: S,
    node :: SupportNode{T},
    stats :: TreeStats;
    fullfilter :: Bool = false,
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false,
    is_singular :: Ref{Bool} = Ref(false)
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    stats.reduction_steps += 1
    for (i, child) in node.children
        if g[i] > 0 || (negative && g[i] < 0) || (fullfilter && g[i] != 0)
            #Look for reducer recursively in the children of this node
            reducer, found_reducer = find_reducer(
                g, gb, child, stats, fullfilter=fullfilter, skipbinomial=skipbinomial,
                negative=negative, is_singular=is_singular
            )
            if found_reducer
                return reducer, found_reducer
            end
        end
    end
    #Check if some binomial in this node works as reducer
    for reducer in node.binomial_list
        # If the reducer is the same element of the GB that was passed as
        #parameter to be skipped, skip it.
        # This is useful in inter-reductions, where the element should not be
        #used to reduce itself.
        if !isnothing(skipbinomial) && reducer == skipbinomial
            continue
        end
        stats.reducers_checked += 1
        if GBElements.reduces(
            g, node.filter, reducer, gb, fullfilter=fullfilter, negative=negative,
            is_singular=is_singular
        )
            return reducer, true
        end
    end
    #No reducer found.
    return g, false
end

end
