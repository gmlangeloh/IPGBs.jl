module ZDDs

using IPGBs.FastBitSets

struct ZDDNode
    index :: Int
    low :: Int
    high :: Int
end

const TOP = ZDDNode(0, 0, 0)
const BOTTOM = ZDDNode(-1, 0, 0)

is_top(node :: ZDDNode) = node.index == 0
is_bottom(node :: ZDDNode) = node.index == -1

struct ZDD
    nodes :: Vector{ZDDNode}
    root :: ZDDNode
end

function ZDD(sets :: Vector{FastBitSets})

end

end