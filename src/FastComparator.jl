"""
Implements fast comparisons of complex objects by building a structure of integers
representing them.
"""
module FastComparator

export Comparator, update!

#Factor by which range is updated when necessary
const RANGE_UPDATE_FACTOR = 2

#Initial integer range for a Comparator
const INITIAL_RANGE = 100000

#How much of the Comparator's range has to be occupied for a rebuild to happen
const OCCUPATION_THRESHOLD = 0.5

"""
A structure allowing for O(1) comparison of elements of `data`. Updating the
structure with new elements of data can take up to O(n). This is done by
associating an integer to each element of `data`, such that data[i] < data[j] iff
magnitudes[i] < magnitudes[j].

TODO Maybe it is possible to optimize insertion to O(log n). Think about it later.
"""
mutable struct Comparator{T} <: Base.Ordering
    data :: Vector{T}
    magnitudes :: Vector{Int}
    sorted_permutation :: Vector{Int}
    range :: Int
    size :: Int

    function Comparator{T}(
        data :: Vector{T}
        range = INITIAL_RANGE
    ) where {T}
        sorted = sortperm(data)
        size = length(data)
        magnitudes = Vector{Int}(undef, size)
        build!(magnitudes, sorted, range, size)
        new{T}(data, magnitudes, sorted, range, size)
    end
end

function Comparator{T}() where {T}
    return Comparator{T}(T[])
end

magnitude(comp :: Comparator{T}, i :: Int) where {T} = comp.magnitudes[i]

"""
Returns the i-th element of the (sorted) data in `comp`.
"""
function ith_element(
    comp :: Comparator{T},
    i :: Int
) :: T where {T}
    return comp.data[comp.sorted_permutation[i]]
end

"""
Ratio of occupied slots in [-comp.range, comp.range]
"""
function occupation_ratio(
    comp :: Comparator{T}
) :: Float64 where {T}
    return comp.size / (2 * comp.range)
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
function build!(
    magnitudes :: Vector{Int},
    sorted :: Vector{Int},
    range :: Int,
    size :: Int
)
    for i in 1:size
        magnitudes[sorted[i]] = magnitude(i, size, range)
    end
end

"""
Returns true iff the element of index `index1` is smaller than that of index
`index2` in comp.data. Uses Comparator logic for comparison in O(1)
"""
function Base.lt(
    comp :: Comparator{T},
    index1 :: Int,
    index2 :: Int
) :: Bool where {T}
    @assert index1 <= comp.size && index2 <= comp.size
    return comp.magnitudes[index1] < comp.magnitudes[index2]
end

"""
Increases range and rebuilds `comp`, equally spacing the magnitudes of the elements
in the new increased range.
"""
function rebuild!(
    comp :: Comparator{T}
) where {T}
    increase_range!(comp)
    build!(comp.magnitudes, comp.sorted_permutation, comp.range, comp.size)
end

"""
Do a binary search for element on the data of `comp` to find the position where
it should be inserted.
"""
function find_position(
    element :: T,
    comp :: Comparator{T}
) :: Int where {T}
    start_index = 1
    end_index = length(comp.data)
    while start_index < end_index
        mid_index = Int(floor((start_index + end_index) / 2))
        pivot = ith_element(comp, mid_index)
        if pivot < element
            start_index = mid_index + 1
        else #element <= pivot
            end_index = mid_index
        end
    end
    position = start_index + 1
    return position
end

"""
Increases the comparator range for magnitudes by RANGE_UPDATE_FACTOR when possible.
If this is not possible, simply increases it to the maximum range possible
(full integers).
"""
function increase_range!(
    comp :: Comparator{T}
) where {T}
    if comp.range == typemax(Int)
        error("Unable to increase FastComparator range anymore.")
    elseif comp.range >= typemax(Int) / RANGE_UPDATE_FACTOR
        comp.range = typemax(Int)
    else
        comp.range *= RANGE_UPDATE_FACTOR
    end
end

"""
Inserts the magnitude of the element of rank `rank` in the magnitude list of
`comp`. May rebuild the whole structure if there are no available magnitudes.
"""
function update_magnitude!(
    comp :: Comparator{T},
    rank :: Int
) where {T}
    #Corner cases: this element is the new smallest / largest element
    if rank == 1
        #New element is the smallest element, thus it is necessary to update the
        #range
        increase_range!(comp)
        mag_new = -comp.range
    elseif rank == comp.size
        #New element is the largest element, and it is necessary to update the
        #range
        increase_range!(comp)
        mag_new = comp.range
    else #Usual case: new element is neither the smallest nor largest
        mag_previous = magnitude(comp, ith_element(comp, rank - 1))
        mag_next = magnitude(comp, ith_element(comp, rank + 1))
        if mag_previous != mag_next
            #New magnitude is the average of its neighbors
            mag_new = Int(round((mag_previous + mag_next) / 2))
        else
            #Problem: the range is too small. It is necessary to rebuild the
            #whole Comparator!
            push!(comp.magnitudes, 0) #Just need to increase the size of the
            #magnitudes vector. We'll get the right magnitude when rebuilding
            rebuild!(comp) #This already sets the magnitude of the new element
            return
        end
    end
    push!(comp.magnitudes, mag_new)
end

"""
Updates `comp` with the magnitudes and sorted positions of any new elements added
to comp.data. May rebuild the whole structure whenever necessary.

TODO is it possible to lower the complexity of this function from O(n) to O(log n)
by avoiding the insertions in a vector? Maybe using some other kind of structure
to store comp.sorted_permutation
"""
function update!(
    comp :: Comparator{T}
) where {T}
    #Easy case: just start to consider the new elements in data
    start = comp.size + 1
    for i in start:length(comp.data)
        position = find_position(comp.data[i], comp)
        insert!(comp.sorted_permutation, position, i)
        comp.size += 1
        update_magnitude!(comp, position)
    end
    @assert comp.size == length(comp.data)
    #For efficiency, rebuild the structure when it starts to get full
    if occupation_ratio(comp) > OCCUPATION_THRESHOLD
        rebuild!(comp)
    end
end

end
