"""
Implements fast comparisons of complex objects by building a structure of integers
representing them.

TODO currently repeated elements give different magnitudes. This is not
necessarily a problem yet, but I should fix that if it ever becomes a
problem
"""
module FastComparator

export Comparator, compare

#Factor by which range is updated when necessary
const RANGE_UPDATE_FACTOR = 2

#How much larger should the range be than the current size of the structure
const NEW_INITIALIZATION_FACTOR = 10000

#How much to increment the range by in corner cases
const RANGE_INCREMENT_SIZE = 50

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
mutable struct Comparator{T, F <: Function} <: Base.Ordering
    data :: Vector{T}
    magnitudes :: Vector{Int}
    sorted_permutation :: Vector{Int}
    range :: Int
    size :: Int
    less_than :: F

    function Comparator{T, F}(
        data :: Vector{T},
        less_than :: F = Base.lt;
        range = INITIAL_RANGE
    ) where {T, F}
        sorted = sortperm(data, lt=less_than)
        size = length(data)
        magnitudes = Vector{Int}(undef, size)
        build!(magnitudes, sorted, range, size)
        new{T, F}(data, magnitudes, sorted, range, size, less_than)
    end
end

function Comparator{T, F}() where {T, F}
    return Comparator{T, F}(T[])
end

magnitude(comp :: Comparator{T, F}, i :: Int) where {T, F} = comp.magnitudes[comp.sorted_permutation[i]]

"""
Returns the i-th element of the (sorted) data in `comp`.
"""
function ith_element(
    comp :: Comparator{T, F},
    i :: Int
) :: T where {T, F}
    return comp.data[comp.sorted_permutation[i]]
end

"""
Ratio of occupied slots in [-comp.range, comp.range]
"""
function occupation_ratio(
    comp :: Comparator{T, F}
) :: Float64 where {T, F}
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
`index2` in comp.data. Uses Comparator logic for comparison in O(1).
"""
function Base.lt(
    comp :: Comparator{T, F},
    index1 :: Int,
    index2 :: Int
) :: Bool where {T, F}
    @assert index1 <= comp.size && index2 <= comp.size
    return comp.magnitudes[index1] < comp.magnitudes[index2]
end

"""
Returns :lt, :eq, :gt if the element indexed by `index1` is smaller, equal or
greater than the one indexed by `index2` respectively. Uses Comparator logic
for comparison in O(1).
"""
function compare(
    comp :: Comparator{T, F},
    index1 :: Int,
    index2 :: Int
) :: Symbol where {T, F}
    c = cmp(comp.magnitudes[index1], comp.magnitudes[index2])
    if c == -1
        return :lt
    elseif c == 1
        return :gt
    end
    return :eq
end

"""
Increases range and rebuilds `comp`, equally spacing the magnitudes of the elements
in the new increased range.
"""
function rebuild!(
    comp :: Comparator{T, F}
) where {T, F}
    increase_range!(comp)
    build!(comp.magnitudes, comp.sorted_permutation, comp.range, comp.size)
end

"""
Do a binary search for element on the data of `comp` to find the position where
it should be inserted.
"""
function find_position(
    element :: T,
    comp :: Comparator{T, F}
) :: Int where {T, F}
    start_index = 1
    end_index = length(comp.data)
    while start_index < end_index
        mid_index = Int(floor((start_index + end_index) / 2))
        pivot = ith_element(comp, mid_index)
        if comp.less_than(pivot, element)
            start_index = mid_index + 1
        else #element <= pivot
            end_index = mid_index
        end
    end
    return start_index
end

"""
Increases the comparator range for magnitudes by RANGE_UPDATE_FACTOR when possible.
If this is not possible, simply increases it to the maximum range possible
(full integers).
"""
function increase_range!(
    comp :: Comparator{T, F};
    use_increment :: Bool = false
) where {T, F}
    if comp.range == typemax(Int)
        if comp.size == typemax(Int)
            error("Unable to increase FastComparator range anymore.")
        end
        #Get a new range that is more adequate to the number of elements in the structure
        if comp.size <= typemax(Int) / NEW_INITIALIZATION_FACTOR
            comp.range = comp.size * NEW_INITIALIZATION_FACTOR
        elseif comp.size < typemax(Int) / 2
            comp.range = comp.size * 2
        end #Otherwise just keep the same range, but this really won't happen with
        #current 64-bit integers as magnitudes
        build!(comp.magnitudes, comp.sorted_permutation, comp.range, comp.size)
    elseif use_increment && comp.range >= typemax(Int) - RANGE_INCREMENT_SIZE
        comp.range = typemax(Int)
    elseif use_increment #and range can be incremented
        comp.range += RANGE_INCREMENT_SIZE
    elseif comp.range >= typemax(Int) / RANGE_UPDATE_FACTOR #!use_increment
        comp.range = typemax(Int)
    else #update by a factor
        comp.range *= RANGE_UPDATE_FACTOR
    end
end

"""
Inserts the magnitude of the element of rank `rank` in the magnitude list of
`comp`. May rebuild the whole structure if there are no available magnitudes.
"""
function update_magnitude!(
    comp :: Comparator{T, F},
    rank :: Int
) where {T, F}
    #Corner cases: this element is the new smallest / largest element
    push!(comp.magnitudes, 0) #Just need to increase the size of the
    #magnitudes vector. We'll set the correct magnitude later in this function
    if rank == 1
        #New element is the smallest element, thus it is necessary to update the
        #range
        increase_range!(comp, use_increment=true)
        mag_new = -comp.range
    elseif rank == comp.size
        #New element is the largest element, and it is necessary to update the
        #range
        increase_range!(comp, use_increment=true)
        mag_new = comp.range
    else #Usual case: new element is neither the smallest nor largest
        mag_previous = magnitude(comp, rank - 1)
        #The next element currently has rank `rank`. It will change after this
        #function terminates
        mag_next = magnitude(comp, rank + 1)
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
    comp.magnitudes[length(comp.magnitudes)] = mag_new
end

"""
Updates `comp` with the magnitudes and sorted positions of any new elements added
to comp.data. May rebuild the whole structure whenever necessary.

TODO is it possible to lower the complexity of this function from O(n) to O(log n)
by avoiding the insertions in a vector? Maybe using some other kind of structure
to store comp.sorted_permutation
"""
function update!(
    comp :: Comparator{T, F}
) where {T, F}
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

"""
Checks whether `comp` satisfies the increasing property, i.e. its elements
in the order given by comp.sorted_permutation are in increasing order

This is useful for testing and debugging.
"""
function check_increasing_property(
    comp :: Comparator{T, F};
    print_data :: Bool = false
) where {T, F}
    if print_data
        println("Checking Comparator property")
        println(ith_element(comp, 1), " ", magnitude(comp, 1))
    end
    for i in 1:(comp.size - 1)
        i_elem = ith_element(comp, i)
        ip1_elem = ith_element(comp, i + 1)
        if print_data
            println(ip1_elem, " ", magnitude(comp, i+1))
        end
        @assert comp.less_than(i_elem, ip1_elem) || isequal(i_elem, ip1_elem)
    end
end

end
