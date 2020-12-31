"""
Implements a combinatorial Buchberger algorithm for Integer Programming
where all data is non-negative. Based on Thomas and Weismantel (1997).
"""
module Buchberger
export buchberger

using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.SupportTrees
using TimerOutputs

"""
Reduces each element of the GB by the previous elements.

TODO count number of removed elements so we can decrement the iteration
index in the main loop
"""
function auto_reduce(
    gb :: Vector{GradedBinomial}
)
    for i in length(gb):-1:1
        nothing
    end
end

"""
Updates gb to a minimal Gröbner Basis.

For more info on minimal GBs, see Lemma 3 from Cox, Little, O'Shea Chapter 2.7.
In summary, it shows that one can remove from a GB any g such that LT(g) is a
multiple of LT(h), for h distinct from g in the GB.
"""
function minimal_basis!(
    gb :: Vector{GradedBinomial},
    tree :: SupportTree{GradedBinomial}
)
    for i in length(gb):-1:1
        g = gb[i]
        if find_reducer(g, gb, tree, skipbinomial=g) != nothing
            deleteat!(gb, i)
            removebinomial!(tree, gb[i])
        end
    end
end

"""
Updates gb to a reduced Gröbner Basis.

TODO bugfix this, it doesn't terminate (at least not in reasonable time)
"""
function reduced_basis!(
    gb :: Vector{GradedBinomial},
    tree :: SupportTree{GradedBinomial}
)
    minimal_basis!(gb, tree)
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        #println("Reducing new element ", i)
        while reducing #&& it < 100
            #TODO I'm getting wrong reducers here
            h = find_reducer(g, gb, tree, negative=true)
            @show g GBElements.fullform(g)
            if h != nothing
                @show h GBElements.fullform(h)
            else
                @show h
            end
            if h != nothing
                GBElements.reduce!(g, h, negative=true) #TODO I think that's it? Check.
                println("new g")
                @show g GBElements.fullform(g)
            else
                reducing = false
            end
        end
        println("post reduction")
        @show g
    end
end

"""
Adds a new element to the current GB.
"""
function update_basis!(
    B :: Vector{GradedBinomial},
    reducer :: SupportTree{GradedBinomial},
    positive_supports :: Vector{FastBitSet},
    negative_supports :: Vector{FastBitSet},
    r :: GradedBinomial
)
    if r.cost < 0
        GBElements.opposite!(r)
    end
    push!(B, r)
    #Update support arrays
    p, n = GBElements.supports(r)
    push!(positive_supports, p)
    push!(negative_supports, n)
    addbinomial!(reducer, B[length(B)])
end

"""
In a maximization problem, if (i, j) ..

TODO actually document this thing, it's not that trivial
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    positive_supports :: Vector{FastBitSet},
    negative_supports :: Vector{FastBitSet}
) :: Bool
    return disjoint(negative_supports[i], negative_supports[j]) ||
        !disjoint(positive_supports[i], positive_supports[j])
end

function supports(
    B :: Vector{GradedBinomial}
) :: Tuple{Vector{FastBitSet}, Vector{FastBitSet}}
    pos_supps = FastBitSet[]
    neg_supps = FastBitSet[]
    for g in B
        p, n = GBElements.supports(g)
        push!(pos_supps, p)
        push!(neg_supps, n)
    end
    return pos_supps, neg_supps
end

"""
Computes a test set / Gröbner Basis for the IP:
max C^T * x
s.t. A * x <= b
x <= u
"""
function buchberger(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    auto_reduce_freq :: Int = 2500
) :: Vector{Vector{Int}}
    n = size(A, 2)
    gb = [ lattice_generator(i, A, C) for i in 1:n ]
    gb = filter(g -> isfeasible(g, A, b, u), gb)
    positive_supports, negative_supports = supports(gb)
    reducer = support_tree(gb, fullfilter=true)
    i = 1
    iteration_count = 0
    spair_count = 0
    zero_reductions = 0
    while i <= length(gb)
        for j in 1:(i-1)
            iteration_count += 1
            if is_support_reducible(i, j, positive_supports, negative_supports)
                continue
            end
            v = gb[i]
            w = gb[j]
            if v.cost < w.cost
                r = w - v #TODO these are relatively expensive
            elseif w.cost < v.cost
                r = v - w
            else #w.cost == v.cost
                if GBElements.lt_tiebreaker(v, w)
                    r = w - v
                else
                    r = v - w
                end
            end
            if isfeasible(r, A, b, u)
                spair_count += 1
                SupportTrees.reduce!(r, gb, reducer)
                if GBElements.iszero(r)
                    zero_reductions += 1
                    continue
                end
                update_basis!(gb, reducer, positive_supports,
                              negative_supports, r)
            end
            if iteration_count % auto_reduce_freq == 0
                auto_reduce(gb)
            end
        end
        i += 1
    end
    @info "Buchberger: S-binomials reduced" iteration_count spair_count zero_reductions
    minimal_basis!(gb, reducer)
    return [ -GBElements.fullform(g) for g in gb ]
end

end
