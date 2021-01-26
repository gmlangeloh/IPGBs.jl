"""
Implements the S-pair triangles described in Roune and Stillman (2012), Section
4.3.
"""
module TriangleHeaps

export TriangleHeap, push_batch!

using DataStructures
using IPGBs.BinomialSets
using IPGBs.SignaturePolynomials

const SigHeap{T} = BinaryHeap{SignaturePair, ModuleMonomialOrdering{T}}

"""
Stores a heap of SignaturePairs (i, j) containing only the S-pair of smallest
signature for each i. Then, for each i, a row of `triangle` contains all the values
of j such that (i, j) is still in the queue, without storing the signatures of
these S-pairs. This reduces memory consumption.

For further improvements in memory consumption, we use the type parameter I to
allow for smaller unsigned integer types. Using Uint16 allows computing with
GBs of up to 65535 elements.
"""
struct TriangleHeap{T, I <: Unsigned}
    heap :: SigHeap{T}
    triangle :: Vector{Vector{I}}
    order :: ModuleMonomialOrdering{T}
    basis :: SigBasis{T}

    function TriangleHeap{T, I}(
        basis :: SigBasis{T},
        order :: ModuleMonomialOrdering{T}
    ) where {T, I}
        heap = BinaryHeap{SignaturePair}(order, [])
        triangle = Vector{I}[]
        new{T, I}(heap, triangle, order, basis)
    end
end

"""
Puts a batch of SignaturePairs (all coming from pairs (i, j) with fixed j) in
`heap`, adding a new row to its triangle. Has to sort the S-pairs by signature,
so its complexity is O(S * n log n) where S is the complexity of signature
comparison and n = length(batch).
"""
function push_batch!(
    heap :: TriangleHeap{T, I},
    batch :: Vector{SignaturePair}
) where {T, I}
    if isempty(batch)
        #Need to add this batch's row to the triangle regardless, to make sure
        #the indices of the next batches work
        push!(heap.triangle, [])
        return
    end
    sort!(batch, order=heap.order) #TODO not sure at all that this works
    #Element of smallest signature goes to the heap...
    push!(heap.heap, batch[1])
    #...the rest go to the triangle, in order of increasing signature
    push!(heap.triangle, Vector{I}(undef, length(batch) - 1))
    j = batch[1].j
    for k in 2:length(batch)
        i = batch[k].i
        heap.triangle[j][k - 1] = i
    end
end

"""
Returns the SignaturePair of lowest signature in `heap`. Assumes that `heap` is
non-empty.
"""
function Base.pop!(
    heap :: TriangleHeap{T, I}
) :: SignaturePair where {T, I}
    #1. Get pair from heap.heap
    pair = pop!(heap.heap)
    #2. Put the new pair from the correct row of heap.triangle in heap.heap
    row = pair.j
    if isempty(heap.triangle[row]) #Nothing to do, no pair available to replace
        return pair
    end
    i = popfirst!(heap.triangle[row])
    new_pair = regular_spair(Int(i), row, heap.basis) #Build the new pair again
    @assert !isnothing(new_pair)
    push!(heap.heap, new_pair)
    return pair
end

"""
Returns true iff `heap` is empty. Just delegates to its internal heap, which is
only empty when `heap` is.
"""
function Base.isempty(
    heap :: TriangleHeap{T, I}
) :: Bool where {T, I}
    return isempty(heap.heap)
end

end
