"""
Reimplementation of 4ti2's FilterReduction tree structure used to store
binomial sets for efficient reductions. Called SupportTrees in reference
to Roune and Stillman (2012).
"""
module SupportTrees
export SupportTree, find_reducer, reduce!, support_tree, addbinomial!,
    removebinomial!

using IPGBs.GBElements

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
    SupportNode{T}() where {T <: AbstractVector} = new(Int[], T[], Int[])
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
    fullfilter :: Bool
    SupportTree{T}(fullfilter) where {T <: AbstractVector{Int}} = begin
        new(SupportNode{T}(), 0, fullfilter)
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
    tree = SupportTree{typeof(gb[1])}(fullfilter)
    for i in 1:length(gb)
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
            newnode = SupportNode{typeof(binomial)}()
            push!(current.children, (i, newnode))
            current = last(current.children)[2]
            #break #I think this is enough, it leaves the loop anyway
        end
    end
    push!(current.binomial_list, binomial)
    if isempty(current.filter) #Set current node's filter if it wasn't already
        for i in binomial_filter #current.filter = f
            push!(current.filter, i)
        end
    end
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
Checks whether `g` is reducible by some element of `gb` and, if so, returns
a reducer. Otherwise, returns `nothing`.
"""
function find_reducer(
    g :: T,
    gb :: S,
    tree :: SupportTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Union{T, Nothing} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    return find_reducer(
        g, gb, tree.root, fullfilter=tree.fullfilter, skipbinomial=skipbinomial,
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
    node :: SupportNode{T};
    fullfilter :: Bool = false,
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Union{T, Nothing} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    for (i, child) in node.children
        if g[i] > 0 || (negative && g[i] < 0) || (fullfilter && g[i] != 0)
            #Look for reducer recursively in the children of this node
            reducer = find_reducer(
                g, gb, child, fullfilter=fullfilter, skipbinomial=skipbinomial,
                negative=negative
            )
            if reducer != nothing #Found a reducer, return it
                return reducer
            end
        end
    end
    #Check if some binomial in this node works as reducer
    for reducer in node.binomial_list
        # If the reducer is the same element of the GB that was passed as
        #parameter to be skipped, skip it.
        # This is useful in inter-reductions, where the element should not be
        #used to reduce itself.
        if skipbinomial != nothing && reducer === skipbinomial
            continue
        end
        if GBElements.reduces(
            g, node.filter, reducer, gb, fullfilter=fullfilter,
            negative=negative)
            return reducer
        end
    end
    return nothing
end

"""
Fully reduce `binomial` by `gb` in-place, finding its normal form. Uses `tree`
to speed up the search for reducers.

`binomial` can also be a monomial.

If `reduction_count` is given, the number of times each reducer was used
in this reduction process is added to reduction_count. This greatly slows
down execution, and is only meant for experimental purposes. The parameter
should be set to `nothing` in practical runs.
"""
function reduce!(
    binomial :: T,
    gb :: S,
    tree :: SupportTree{T};
    reduction_count :: Union{Vector{Int}, Nothing} = nothing,
    skipbinomial :: Union{T, Nothing} = nothing
) :: Bool where {T <: GBElement, S <: AbstractVector{T}}
    reducer = find_reducer(binomial, gb, tree, skipbinomial=skipbinomial)
    reduced_to_zero = false
    while reducer != nothing
        reduced_to_zero = GBElements.reduce!(binomial, reducer)
        if reduction_count != nothing
            for i in 1:length(gb)
                if gb[i] === reducer
                    reduction_count[i] += 1
                    break
                end
            end
        end
        if reduced_to_zero
            break
        end
        reducer = find_reducer(binomial, gb, tree, skipbinomial=skipbinomial)
    end
    return reduced_to_zero
end

end
