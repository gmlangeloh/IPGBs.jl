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
Returns the i-th element of the (sorted) data in `comp`.

TODO should I make Comparator implement the array interface?
Then I could just get the element with indexing
"""
function ith_element(
    comp :: Comparator{T, O},
    i :: Int
) :: T where {T, O <: Base.Ordering}
    return comp.data[comp.sorted_permutation[i]]
end

"""
Ratio of occupied slots in [-comp.range, comp.range]
"""
function occupation_ratio(
    comp :: Comparator{T, O}
) :: Float64 where {T, O <: Base.Ordering}
    return comp.size / (2 * comp.range)
end

function magnitude(
    comp :: Comparator{T, O},
    index :: Int
) :: Int where {T, O <: Base.Ordering}
    return comp.magnitudes[index]
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
Inserts the magnitude of the element of rank `rank` in the magnitude list of
`comp`. May rebuild the whole structure if there are no available magnitudes.
"""
function update_magnitude!(
    comp :: Comparator{T, O},
    rank :: Int
) where {T, O <: Base.Ordering}
    #Corner cases: this element is the new smallest / largest element
    if rank == 1
        #New element is the smallest element, thus it is necessary to update the
        #range

        #TODO update range and compute new magnitude
    elseif rank == comp.size
        #New element is the largest element, and it is necessary to update the
        #range

        #TODO update range and compute new magnitude
    else #Usual case: new element is neither the smallest nor largest
        mag_previous = magnitude(comp, ith_element(comp, rank - 1))
        mag_next = magnitude(comp, ith_element(comp, rank + 1))
        if mag_previous != mag_next
            #New magnitude is the average of its neighbors
            mag_new = Int(round((mag_previous + mag_next) / 2))
        else
            #Problem: the range is too small. It is necessary to rebuild the
            #whole Comparator!
            rebuild!(comp) #This already sets the magnitude of the new element
            return
        end
    end
    push!(comp.magnitudes, mag_new)
end

"""
Updates `comp` with the magnitudes and sorted positions of any new elements added
to comp.data. May rebuild the whole structure whenever necessary.
"""
function update!(
    comp :: Comparator{T, O}
) where {T, O <: Base.Ordering}
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
