"""
Implements a combinatorial Buchberger algorithm for Integer Programming
where all data is non-negative. Based on Thomas and Weismantel (1997).
"""
module Buchberger
export buchberger

using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.SupportTrees

using TimerOutputs

"""
Reduces each element of the GB by the previous elements.

TODO count number of removed elements so we can decrement the iteration
index in the main loop
"""
function auto_reduce(
    gb :: Vector{T}
) where {T <: GBElement}
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
    gb :: Vector{T},
    tree :: SupportTree{T}
) where  {T <: GBElement}
    for i in length(gb):-1:1
        g = gb[i]
        red = find_reducer(g, gb, tree, skipbinomial=g)
        if red != nothing
            @show red
            @show g
            deleteat!(gb, i)
            removebinomial!(tree, g)
        end
    end
end

"""
Updates gb to a reduced Gröbner Basis.

TODO bugfix this, it doesn't terminate or is just buggy overall
"""
function reduced_basis!(
    gb :: Vector{T},
    tree :: SupportTree{T}
) where {T <: GBElement}
    minimal_basis!(gb, tree)
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        while reducing
            h = find_reducer(g, gb, tree, negative=true)
            if h != nothing
                GBElements.reduce!(g, h, negative=true)
            else
                reducing = false
            end
        end
    end
end

"""
Adds a new element to the current GB.
"""
function update_basis!(
    B :: Vector{T},
    reducer :: SupportTree{T},
    positive_supports :: Vector{FastBitSet},
    negative_supports :: Vector{FastBitSet},
    r :: T
) where {T <: GBElement}
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
Returns true if (i, j) should be discarded.

In a maximization problem, if (i, j) ..

TODO actually document this thing, it's not that trivial
"""
function is_support_reducible(
    i :: Int,
    j :: Int,
    positive_supports :: Vector{FastBitSet},
    negative_supports :: Vector{FastBitSet},
    minimization :: Bool
) :: Bool
    if minimization
        return !disjoint(negative_supports[i], negative_supports[j]) ||
            disjoint(positive_supports[i], positive_supports[j])
    end
    #Maximization problem
    return disjoint(negative_supports[i], negative_supports[j]) ||
        !disjoint(positive_supports[i], positive_supports[j])
end

function supports(
    B :: Vector{T}
) :: Tuple{Vector{FastBitSet}, Vector{FastBitSet}} where {T <: GBElement}
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
Builds the S-binomial given by gb[i] and gb[j].
"""
function build_sbin(
    i :: Int,
    j :: Int,
    gb :: Vector{T}
) :: T where {T <: AbstractVector{Int}}
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
    return r
end

"""
Computes a test set / Gröbner Basis for the IP:
max C^T * x
s.t. A * x <= b
x <= u

Structure refers to the data structure used to represent binomials internally.
It can be either `Binomial` or `GradedBinomial`.
"""
function buchberger(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    structure :: DataType = Binomial,
    auto_reduce_freq :: Int = 2500
) :: Vector{Vector{Int}}
    @assert structure == Binomial || structure == GradedBinomial
    n = size(A, 2)
    if structure == Binomial
        #The reductions without fullfilter only work correctly if the problem
        #is in minimization form. Thus we take the opposite of C instead, as
        #this is easier than changing everything else
        C = -C
        minimization = true
        gb = [ lattice_generator_binomial(i, A, C) for i in 1:n ]
        A, b, C, u = GBTools.normalize(A, b, C, u)
    else
        minimization = false
        gb = [ lattice_generator_graded(i, A, C) for i in 1:n ]
    end
    gb = Base.filter(g -> isfeasible(g, A, b, u), gb)
    positive_supports, negative_supports = supports(gb)
    reducer = support_tree(gb, fullfilter=(structure == GradedBinomial))
    i = 1
    iteration_count = 0
    spair_count = 0
    zero_reductions = 0
    #Main loop: generate all relevant S-binomials and reduce them
    while i <= length(gb)
        for j in 1:(i-1)
            iteration_count += 1
            if is_support_reducible(
                i, j, positive_supports, negative_supports, minimization
            )
                continue
            end
            r = build_sbin(i, j, gb)
            if i == 6 && j == 1
                @show r
            end
            if isfeasible(r, A, b, u)
                spair_count += 1
                reduced_to_zero = SupportTrees.reduce!(r, gb, reducer)
                if reduced_to_zero
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
    #Convert the basis to the same format 4ti2 uses
    if structure == Binomial
        minimal_basis!(gb, reducer)
        output_basis = gb
    else
        minimal_basis!(gb, reducer)
        output_basis = [ -GradedBinomials.fullform(g) for g in gb ]
    end
    return output_basis
end

end
