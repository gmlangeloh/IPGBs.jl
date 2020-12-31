"""
Specialized module for faster bitsets in the context of GBs / IP.
"""
module FastBitSets
export FastBitSet, iszero, intersect, disjoint, makebitset

mutable struct FastBitSet{T <: Integer} <: AbstractVector{Bool}
    content :: T
    length :: Int
end

function makebitset(
    length :: Int
)
    if length <= 64
        return FastBitSet{Int64}(0, 64)
    elseif length <= 128
        return FastBitSet{Int128}(0, 128)
    end
    return FastBitSet{BigInt}(0, length)
end

function makebitset(
    length :: Int,
    indices :: Array{Int}
)
    bitset = makebitset(length)
    for i in indices
        bitset[i] = true
    end
    return bitset
end

#
# Implementation of the AbstractVector interface
#

function Base.size(
    bitset :: FastBitSet
) :: Tuple
    return (bitset.length, )
end

function Base.getindex(
    bitset :: FastBitSet,
    i :: Int
) :: Bool
    mask = 1 << (i - 1)
    return mask & bitset.content == 0 ? false : true
end

function Base.setindex!(
    bitset :: FastBitSet,
    v :: Bool,
    i :: Int
)
    if v
        if !bitset[i]
            bitset.content += 1 << (i - 1)
        end
    else
        if bitset[i]
            bitset.content -= 1 << (i - 1)
        end
    end
end

function Base.length(
    bitset :: FastBitSet
) :: Int
    return bitset.length
end

#
# Additional operations on bitsets.
#

function intersect(
    bitset1 :: FastBitSet,
    bitset2 :: FastBitSet
) :: FastBitSet
    return FastBitSet(bitset1.content & bitset2.content, bitset1.length)
end

function isempty(
    bitset :: FastBitSet
) :: Bool
    return bitset.content == 0
end

function disjoint(
    bitset1 :: FastBitSet,
    bitset2 :: FastBitSet
) :: Bool
    return isempty(intersect(bitset1, bitset2))
end

end
