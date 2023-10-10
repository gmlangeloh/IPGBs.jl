"""
Reimplementation of 4ti2's FilterReduction tree structure used to store
binomial sets for efficient reductions. Called SupportTrees in reference
to Roune and Stillman (2012).
"""
module SupportTrees
export ReductionTree, SupportTree, CacheTree, find_reducer, support_tree, add_binomial!,
    remove_binomial!, enumerate_reducers, find_reducer_iter


using LRUCache

using IPGBs
using IPGBs.GBElements
using IPGBs.Statistics

abstract type ReductionTree{T} end

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
mutable struct SupportTree{T <: AbstractVector{Int}} <: ReductionTree{T}
    root :: SupportNode{T}
    size :: Int
    depth :: Int #This is useful to know how balanced the tree is
    fullfilter :: Bool
    stats :: TreeStats

    SupportTree{T}(fullfilter) where {T <: AbstractVector{Int}} = begin
        new(SupportNode{T}(), 0, 1, fullfilter, TreeStats())
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
function SupportTree{T}(
    gb :: S;
    fullfilter :: Bool = false
) where { T <: AbstractVector{Int}, S <: AbstractVector{T}}
    tree = SupportTree{T}(fullfilter)
    for g in gb
        add_binomial!(tree, g)
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
function add_binomial!(
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
    tree.size += 1
end

"""
Removes `binomial` from `tree` by doing a lookup operation and then just
removes it from the corresponding binomial_list.

If `binomial` is not in `tree`, nothing will happen.
"""
function remove_binomial!(
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
        if GBElements.reduces(g, node.filter, reducer, gb, negative=negative)
            push!(reducers, reducer)
        end
    end
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
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    #The fullfilter parameter was removed as it affects performance and
    #it was only useful when T == GradedBinomial.
    return find_reducer(
        g, gb, tree.root, tree.stats, skipbinomial=skipbinomial, 
        negative=negative
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
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    stats.reduction_steps += 1
    sign = negative ? -1 : 1
    for (i, child) in node.children
        if sign * g[i] > 0
            #Look for reducer recursively in the children of this node
            reducer, found_reducer = find_reducer(
                g, gb, child, stats, skipbinomial=skipbinomial,
                negative=negative
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
            g, node.filter, reducer, gb, negative=negative
        )
            return reducer, true
        end
    end
    #No reducer found.
    return g, false
end

#TODO: Fix the performance of these cache trees. The issue is likely that
#removing an element from the cache_tree keeps many useless internal nodes.
#As most elements are used for reductions at least once, in time the cache
#tree will have the same structure as the full tree, only with fewer
#elements. Thus, searching it will not be efficient, and at that point
#it's better to just ignore it.
# A possible fix is to rebuild this cache tree periodically, or to adjust
#remove_binomial! to remove useless internal nodes.
mutable struct CacheTree{T <: AbstractVector{Int}} <: ReductionTree{T}
    cache_tree :: SupportTree{T}
    full_tree :: SupportTree{T}
    lru :: LRU{T, Bool}
    max_cache_size :: Int
    cache_misses :: Int
    last_evicted :: Vector{T}

    function CacheTree{T}(fullfilter :: Bool = false) where {T <: AbstractVector{Int}}
        cache = SupportTree{T}(fullfilter)
        tree = SupportTree{T}(fullfilter)
        max_cache_size = IPGBs.CACHE_TREE_SIZE
        last_evicted = T[]
        recover_binomial(binomial, val) = push!(last_evicted, binomial)
        lru = LRU{T, Bool}(maxsize=max_cache_size, finalizer=recover_binomial)
        return new(cache, tree, lru, max_cache_size, 0, last_evicted)
    end
end

function CacheTree{T}(
    gb :: S; 
    fullfilter :: Bool = false
) where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    tree = CacheTree{T}(fullfilter)
    for g in gb
        add_binomial!(tree, g)
    end
    return tree
end

function add_binomial!(
    tree :: CacheTree{T}, 
    binomial :: T
) where {T <: AbstractVector{Int}}
    add_binomial!(tree.full_tree, binomial)
    if length(tree.lru) < tree.max_cache_size
        tree.lru[binomial] = true
        add_binomial!(tree.cache_tree, binomial)
    end
end

function remove_binomial!(
    tree :: CacheTree{T}, 
    binomial :: T
) where {T <: AbstractVector{Int}}
    if haskey(tree.lru, binomial)
        delete!(tree.lru, binomial)
        remove_binomial!(tree.cache_tree, binomial)
    end
    remove_binomial!(tree.full_tree, binomial)
end

function update_cache!(
    tree :: CacheTree{T}, 
    binomial :: T
) where {T <: AbstractVector{Int}}
    if !haskey(tree.lru, binomial)
        tree.cache_misses += 1
        tree.lru[binomial] = true
        #If the cache is full, adding the new element to lru will evict
        #a previous element and put it in last_removed. This element has
        #to be removed from the cache_tree as well.
        if !isempty(tree.last_evicted)
            evicted = pop!(tree.last_evicted)
            remove_binomial!(tree.cache_tree, evicted)
        end
        add_binomial!(tree.cache_tree, binomial)
    end
end

function find_reducer(
    g :: T, 
    gb :: S, 
    tree :: CacheTree{T}; 
    skipbinomial :: Union{T, Nothing} = nothing, 
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    #First search for some reducer for binomial in the cache. If it isn't 
    #found there, look for it in the full tree. If it's found in the full 
    #tree, that's a cache miss, so update the cache with the reducer.
    reducer, found_cache_reducer = find_reducer(
        g, gb, tree.cache_tree, skipbinomial=skipbinomial, negative=negative
    )
    if found_cache_reducer
        #Update the number of times the reducer was hit in the cache.
        _ = get(tree.lru, reducer, false)
        return reducer, true
    end
    reducer, found_reducer = find_reducer(
        g, gb, tree.full_tree, skipbinomial=skipbinomial, negative=negative
    )
    if found_reducer
        update_cache!(tree, reducer)
        return reducer, true
    end
    return g, false
end

end
