using IPGBs.FastBitSets
using IPGBs.ZDDs

b1 = FastBitSet(6, [1, 2])
b2 = FastBitSet(6, [2, 4])
b3 = FastBitSet(6, [3, 4])
b4 = FastBitSet(6, [1, 5, 6])
b5 = FastBitSet(6, [1, 3, 5])
b6 = FastBitSet(6, [4, 6])
bs = [b1, b2, b3, b4, b5, b6]

zdd = ZDD(bs)
println(zdd)