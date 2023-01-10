module SupportMatrixTrees
export SupportMatrixTree, addbinomial!, removebinomial!, find_reducer, find_reducer_rec

using ElasticArrays

using IPGBs.GBElements

mutable struct SupportMatrixTree{T <: AbstractVector{Int}}
    matrix_tree :: ElasticMatrix{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}
    max_children :: Int
    binomial_list :: Vector{Vector{T}}
    filter_list :: Vector{Vector{Int}}
    node_stack :: Vector{Tuple{Int, Int}}
    stack_size :: Int
    stats :: String #TODO: Collect better statistics later!!!

    SupportMatrixTree{T}(n :: Int) where {T <: AbstractVector{Int}} = begin
        mat = ElasticMatrix{Tuple{Int, Int}}(undef, n, 0)
        bins = Vector{T}[]
        filters = Vector{Int}[]
        stack = [(0, 0) for _ in 1:n]
        #Keep the same preallocated stack vector for efficiency
        stack_size = 0
        new(mat, n, bins, filters, stack, stack_size, "")
    end 
end

function stack_push!(
    tree :: SupportMatrixTree{T}, 
    elem :: Tuple{Int, Int}
) where {T <: AbstractVector{Int}}
    tree.stack_size += 1
    tree.node_stack[tree.stack_size] = elem
end

function stack_pop!(
    tree :: SupportMatrixTree{T}
) where {T <: AbstractVector{Int}}
    result = tree.node_stack[tree.stack_size]
    tree.stack_size -= 1
    return result
end

stack_isempty(
    tree :: SupportMatrixTree{T}
) where {T <: AbstractVector{Int}} = tree.stack_size == 0

stack_empty!(
    tree :: SupportMatrixTree{T}
) where {T <: AbstractVector{Int}} = tree.stack_size = 0

function add_empty_node!(
    tree :: SupportMatrixTree{T}
) where {T <: AbstractVector{Int}}
    new_node = Vector{Tuple{Int, Int}}()
    sizehint!(new_node, tree.max_children)
    for _ in 1:tree.max_children
        push!(new_node, (0, 0))
    end
    append!(tree.matrix_tree, new_node)
    push!(tree.filter_list, Int[])
    push!(tree.binomial_list, T[])
end

function SupportMatrixTree(
    gb :: S
) :: SupportMatrixTree{T} where { T <: AbstractVector{Int}, S <: AbstractVector{T}}
    @assert !isempty(gb)
    tree = SupportMatrixTree{T}(length(gb[1]))
    add_empty_node!(tree)
    for g in gb
        addbinomial!(tree, g)
    end
    return tree
end

function addbinomial!(
    tree :: SupportMatrixTree{T}, 
    binomial :: T
) where {T <: AbstractVector{Int}}
    current_index = 1 #Column 1 represents the root of the tree
    binomial_filter = GBElements.filter(binomial)
    for i in binomial_filter
        j = 1
        found_i = false
        while j <= tree.max_children
            if tree.matrix_tree[j, current_index][1] == 0
                break
            elseif tree.matrix_tree[j, current_index][1] == i
                found_i = true
                break
            end
            j += 1
        end
        if found_i #Next node in the path already exists
            current_index = tree.matrix_tree[j, current_index][2]
        else #Next node doesn't exist, create it
            add_empty_node!(tree)
            new_node_index = size(tree.matrix_tree, 2)
            tree.matrix_tree[j, current_index] = (i, new_node_index)
            current_index = new_node_index
            break
        end
    end
    push!(tree.binomial_list[current_index], binomial)
    if isempty(tree.filter_list[current_index])
        tree.filter_list[current_index] = copy(binomial_filter)
    end
end

function removebinomial!(
    tree :: SupportMatrixTree, 
    binomial :: T
) where {T <: AbstractVector{Int}}
    current_index = 1
    binomial_filter = GBElements.filter(binomial)
    for i in binomial_filter
        j = 1
        found_i = false
        while j <= tree.max_children
            if tree.matrix_tree[j, current_index][1] == 0
                break
            elseif tree.matrix_tree[j, current_index][1] == i
                found_i = true
                break
            end
            j += 1
        end
        if found_i #Next node in the path already exists
            current_index = tree.matrix_tree[j, current_index][2]
        end
    end
    for i in eachindex(tree.binomial_list[current_index])
        g = tree.binomial_list[current_index][i]
        if g == binomial
            deleteat!(tree.binomial_list[current_index], i)
            return
        end
    end
end

function check_node_for_reducers(
    current_index :: Int,
    tree :: SupportMatrixTree{T},
    g :: T,
    gb :: S;
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    for reducer in tree.binomial_list[current_index]
        if !isnothing(skipbinomial) && reducer == skipbinomial
            continue
        end
        f = tree.filter_list[current_index]
        if GBElements.reduces(
            g, f, reducer, gb, negative=negative
        )
            return reducer, true
        end
    end
    return g, false
end

function find_reducer(
    g :: T,
    gb :: S,
    tree :: SupportMatrixTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    sign = negative ? -1 : 1
    reducer = g
    found_reducer = false
    @assert stack_isempty(tree)
    stack_push!(tree, (1, 1))
    while !stack_isempty(tree)
        current_index, i = stack_pop!(tree)
        while i <= tree.max_children
            j, child_index = tree.matrix_tree[i, current_index]
            if j == 0 #This node has no more children, ignore it
                break
            elseif sign * g[j] > 0
                #This node has a child that may contain some reducer,
                #so search it. Store the current node for backtracking
                stack_push!(tree, (current_index, i + 1))
                stack_push!(tree, (child_index, 1))
                break
            end
            i += 1
        end
        reducer, found_reducer = check_node_for_reducers(
            current_index, tree, g, gb, skipbinomial=skipbinomial,
            negative=negative
        )
        if found_reducer
            break
        end
    end
    stack_empty!(tree)
    return reducer, found_reducer
end

#There is no apparent performance advantage in the recursive version nor in the
#iterative one. However, it's easier to profile the iterative version.
function find_reducer_rec(
    g :: T,
    gb :: S,
    tree :: SupportMatrixTree{T};
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    return find_reducer_rec(
        g, gb, tree, 1, skipbinomial=skipbinomial, negative=negative
    ) 
end

function find_reducer_rec(
    g :: T,
    gb :: S,
    tree :: SupportMatrixTree{T},
    node_index :: Int;
    skipbinomial :: Union{T, Nothing} = nothing,
    negative :: Bool = false
) :: Tuple{T, Bool} where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    sign = negative ? -1 : 1
    i = 1
    while i <= tree.max_children
        j, child_index = tree.matrix_tree[i, node_index]
        if j == 0
            break
        end
        if sign * g[j] > 0
            reducer, found_reducer = find_reducer_rec(
                g, gb, tree, child_index, skipbinomial=skipbinomial, negative=negative
            )
            if found_reducer
                return reducer, found_reducer
            end
        end
        i += 1
    end
    for reducer in tree.binomial_list[node_index]
        if !isnothing(skipbinomial) && reducer == skipbinomial
            continue
        end
        if GBElements.reduces(
            g, tree.filter_list[node_index], reducer, gb, negative=negative
        )
            return reducer, true
        end
    end
    return g, false
end

end