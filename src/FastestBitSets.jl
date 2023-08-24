
module FastestBitSets
#Currently slower than FastBitSets

const INT_TYPE = UInt64
const BITS_PER_WORD = 8 * sizeof(INT_TYPE)

struct FastestBitSet
    memory :: Ptr{INT_TYPE}
    size :: Int
end

Base.show(io :: IO, bs :: FastestBitSet) = print(io, bs.memory)
Base.size(bs :: FastestBitSet) = (bs.size,)

word_and_index(i :: Int) = ceil(Int, i / BITS_PER_WORD), i % BITS_PER_WORD

function FastestBitSet(
    vars :: Int,
    indices :: Vector{Int}
)
    words = Int(ceil(vars / BITS_PER_WORD))
    memory = Ptr{INT_TYPE}(Libc.calloc(words, sizeof(INT_TYPE)))
    for i in indices
        word, index = word_and_index(i)
        prev = unsafe_load(memory, word)
        unsafe_store!(memory, prev | (1 << index), word)
    end
    return FastestBitSet(Ptr{INT_TYPE}(memory), words)
end

function disjoint(bs1 :: Ptr{INT_TYPE}, bs2 :: Ptr{INT_TYPE}, size :: Int) :: Bool
    for i in 1:size
        if unsafe_load(bs1, i) & unsafe_load(bs2, i) != 0
            return false
        end
    end
    return true
end

disjoint_c = @cfunction(disjoint, Bool, (Ptr{INT_TYPE}, Ptr{INT_TYPE}, Int))

disjoint(bs1 :: FastestBitSet, bs2 :: FastestBitSet) = begin
    @ccall $disjoint_c(bs1.memory :: Ptr{INT_TYPE}, bs2.memory :: Ptr{INT_TYPE}, bs1.size :: Int)::Bool
end

end