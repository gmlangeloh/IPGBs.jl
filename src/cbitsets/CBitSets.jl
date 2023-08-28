module CBitSets

const CLIB = "/home/gmlangeloh/.julia/dev/IPGBs/src/cbitsets/cbitsets.so"

struct CBitSet
    words :: Ptr{UInt}
    size :: UInt
end

function CBitSet(vars :: Int, indices :: Vector{Int})
    bs = @ccall CLIB.bitset_create(vars :: Cuint) :: Ptr{CBitSet}
    @ccall CLIB.bitset_fill(bs :: Ptr{CBitSet}, indices :: Ptr{Cint}, length(indices) :: Cuint) :: Cvoid
    return unsafe_load(bs)
end

function disjoint(
    bs1 :: CBitSet,
    bs2 :: CBitSet
)
    ptr_bs1 = Ref(bs1)
    ptr_bs2 = Ref(bs2)
    return @ccall CLIB.is_disjoint(ptr_bs1 :: Ref{CBitSet}, ptr_bs2 :: Ref{CBitSet}) :: Bool
end

end