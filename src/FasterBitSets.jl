module FasterBitSets

export FasterBitSet, disjoint

const BITS_PER_WORD = 8 * sizeof(Int)

struct FasterBitSet
    data :: Vector{Int}
    words :: Int

    function FasterBitSet(
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

FasterBitSet(vars :: Int) = FasterBitSet(vars, Int[])

word_and_index(i :: Int) = Int(ceil(i / BITS_PER_WORD)), i % BITS_PER_WORD

#
# Implementation of the AbstractVector interface
#

function Base.size(
    bitset :: FasterBitSet
) :: Tuple
    return (bitset.words, )
end

"""
Truth value associated by this bitset to the i-th variable.
"""
function Base.getindex(
    bitset :: FasterBitSet,
    i :: Int
) :: Bool
    word, index = word_and_index(i)
    mask = 1 << index
    return mask & bitset.data[word] == 0 ? false : true
end

function Base.setindex!(
    bitset :: FasterBitSet,
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
    bitset :: FasterBitSet
) :: Int
    return bitset.words
end

#
# Additional operations on bitsets
#

function isempty(
    bitset :: FasterBitSet
) :: Bool
    for i in 1:length(bitset)
        if bitset.data[i] != 0
            return false
        end
    end
    return true
end

function disjoint(
    bitset1 :: FasterBitSet,
    bitset2 :: FasterBitSet
) :: Bool
    for i in 1:length(bitset1)
        if bitset1.data[i] & bitset2.data[i] != 0
            return false
        end
    end
    return true
end

end
