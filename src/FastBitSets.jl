"""
Specialized module for faster bitsets in the context of GBs / IP.
"""
module FastBitSets
export FastBitSet, disjoint, BitTriangle, add_row!

const BITS_PER_WORD = 8 * sizeof(Int)

#
# Fast (uni-dimensional) bitsets that are useful mainly when we know their
# maximum size (e.g. number of variables of an instance).
#

struct FastBitSet
    data :: Vector{Int}
    words :: Int

    function FastBitSet(
        vars :: Int,
        indices :: Vector{Int}
    )
        words = Int(ceil(vars / BITS_PER_WORD))
        data = zeros(Int, words)
        for i in indices
            word, index = word_and_index(i)
            data[word] += 1 << index
        end
        new(data, words)
    end
end

function Base.show(
    io :: IO,
    bitset :: FastBitSet
)
    for i in bitset.words:-1:1
        word = bitset.data[i]
        for j in BITS_PER_WORD:-1:1
            mask = 1 << (j - 1)
            bit = word & mask != 0 ? "1" : "0"
            print(io, bit)
        end
    end
end

FastBitSet(vars :: Int) = FastBitSet(vars, Int[])

word_and_index(i :: Int) = Int(ceil(i / BITS_PER_WORD)), i % BITS_PER_WORD

#
# Implementation of the AbstractVector interface
#

function Base.size(
    bitset :: FastBitSet
) :: Tuple
    return (bitset.words, )
end

"""
Truth value associated by this bitset to the i-th variable.
"""
function Base.getindex(
    bitset :: FastBitSet,
    i :: Int
) :: Bool
    word, index = word_and_index(i)
    mask = 1 << index
    return mask & bitset.data[word] == 0 ? false : true
end

function Base.setindex!(
    bitset :: FastBitSet,
    v :: Bool,
    i :: Int
)
    word, index = word_and_index(i)
    if v
        if !bitset[i]
            bitset.data[word] += 1 << index
        end
    else
        if bitset[i]
            bitset.data[word] -= 1 << index
        end
    end
end

function Base.length(
    bitset :: FastBitSet
) :: Int
    return bitset.words
end

#
# Additional operations on bitsets
#

function isempty(
    bitset :: FastBitSet
) :: Bool
    for i in 1:length(bitset)
        if bitset.data[i] != 0
            return false
        end
    end
    return true
end

function disjoint(
    bitset1 :: FastBitSet,
    bitset2 :: FastBitSet
) :: Bool
    for i in 1:length(bitset1)
        if bitset1.data[i] & bitset2.data[i] != 0
            return false
        end
    end
    return true
end
#
# Bit triangles offering O(1) access to binary data relative to a pair (i, j)
#

"""
A BitTriangle stores a bit of data for each pair (i, j) with i != j. Data can
be set and accessed for both (i, j) and (j, i), even though it is only stored
once per pair.
"""
struct BitTriangle
    data :: Vector{BitVector}

    BitTriangle() = new(BitVector[])
end

"""
Adds a new row (initialized to false) to this BitTriangle. If this is the n-th
row, it will have (n-1) elements.
"""
function add_row!(
    triangle :: BitTriangle
)
    new_row = zeros(Int, length(triangle.data))
    push!(triangle.data, new_row)
end

"""
Checks whether the indices i, j can be a pair in a bit triangle. In practice,
this means they have to be different.
"""
@inline check_indices(i, j) = @assert i != j

"""
Swaps indices so that they are access a valid pair in a bit triangle.
"""
@inline function triangle_indices(i, j)
    if i < j
        return j, i
    end
    return i, j
end

function Base.getindex(
    triangle :: BitTriangle,
    i :: Int,
    j :: Int
) :: Bool
    check_indices(i, j)
    i, j = triangle_indices(i, j)
    return triangle.data[i][j]
end

function Base.setindex!(
    triangle :: BitTriangle,
    value :: Bool,
    i :: Int,
    j :: Int
)
    check_indices(i, j)
    i, j = triangle_indices(i, j)
    triangle.data[i][j] = value
end

function Base.show(
    io :: IO,
    triangle :: BitTriangle
)
    n = length(triangle.data)
    for i in 1:n
        if i < n
            println(io, triangle.data[i])
        else
            print(io, triangle.data[i])
        end
    end
end

end
