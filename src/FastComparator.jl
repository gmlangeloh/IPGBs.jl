"""
Implements fast comparisons of complex objects by building a structure of integers
representing them.
"""
module FastComparator

#Factor by which range is updated when necessary
const RANGE_UPDATE_FACTOR = 2

#Initial integer range for a Comparator
const INITIAL_RANGE = 100000

#How much of the Comparator's range has to be occupied for a rebuild to happen
const OCCUPATION_THRESHOLD = 0.5

#Compares stuff of type T?
mutable struct Comparator{T, O <: Base.Ordering} <: Base.Ordering
    data :: Vector{T}
    order :: O #TODO Maybe I don't need this?
    magnitudes :: Vector{Int}
    sorted_permutation :: Vector{Int}
    range :: Int
    size :: Int

    function Comparator{T, O}(data :: Vector{T}, order :: O,
                              range = INITIAL_RANGE) where {T, O <: Base.Ordering}
        sorted, magnitudes = build(data, range)
        size = length(data)
        new{T, O}(data, order, magnitudes, sorted, range, size)
    end
end

function Comparator{T, O}(order :: O) where {T, O <: Base.Ordering}
    return Comparator{T, O}(T[], order)
end

"""
Computes the magnitude of the i-th smallest element out of n elements with the
given range. The elements are equally spaced in the range.
"""
function magnitude(
    i :: Int,
    n :: Int,
    range :: Int
) :: Int
    step = Int(floor(2 * range / n))
    return -range + i * step
end

"""
Returns a vector of sorted permutation indices and integer magnitudes of elements
of `data` in the range given by [-range, range].
"""
function build(
    data :: Vector{T},
    range :: Int
) :: Tuple{Vector{Int}, Vector{Int}} where {T}
    n = length(data)
    magnitudes = Vector{Int}(undef, n)
    sorted = sortperm(data)
    for i in 1:n
        magnitudes[sorted[i]] = magnitude(i, n, range)
    end
    return sorted, magnitudes
end

function Base.lt(
    comp :: Comparator{T, O},
    index1 :: Int,
    index2 :: Int
) :: Bool where {T, O <: Base.Ordering}
    @assert index1 <= comp.size && index2 <= comp.size
    return comp.magnitudes[index1] < comp.magnitudes[index2]
end

function rebuild!(
    comp :: Comparator{T, O}
) where {T, O <: Base.Ordering}

end

"""
Do a binary search for element on the data of `comp` to find the position where
it should be inserted.
"""
function find_position(
    element :: T,
    comp :: Comparator{T, O}
) :: Int where {T, O <: Base.Ordering}

end

function update!(
    comp :: Comparator{T, O}
) where {T, O <: Base.Ordering}
    #Easy case: just start to consider the new elements in data
    start = comp.size + 1
    for i in start:length(comp.data)

    end

    #Hard case: Has to call rebuild from time to time
end

end
