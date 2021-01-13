"""
Implements a combinatorial Buchberger algorithm for Integer Programming
where all data is non-negative. Based on Thomas and Weismantel (1997).
"""
module Buchberger
export buchberger

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.SupportTrees

"""
Reduces each element of the GB by the previous elements.

TODO count number of removed elements so we can decrement the iteration
index in the main loop
"""
function auto_reduce(
    gb :: BinomialSet{T, S}
) where {T <: GBElement, S <: GBOrder}
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
    gb :: BinomialSet{T, S}
) where  {T <: GBElement, S <: GBOrder}
    for i in length(gb):-1:1
        g = gb[i]
        reducer = find_reducer(g, gb, reduction_tree(gb), skipbinomial=g)
        if !isnothing(reducer)
            @show reducer
            @show g
            deleteat!(gb, i)
            removebinomial!(reduction_tree(gb), g)
        end
    end
end

"""
Updates gb to a reduced Gröbner Basis.

TODO bugfix this, it doesn't terminate or is just buggy overall
"""
function reduced_basis!(
    gb :: BinomialSet{T, S}
) where {T <: GBElement, S <: GBOrder}
    minimal_basis!(gb)
    for i in length(gb):-1:1
        g = gb[i]
        reducing = true
        while reducing
            h = find_reducer(g, gb, reduction_tree(gb), negative=true)
            if !isnothing(h)
                GBElements.reduce!(g, h, negative=true)
            else
                reducing = false
            end
        end
    end
end

"""
Builds the S-binomial given by gb[i] and gb[j].

TODO probably should use MonomialOrder to orientate stuff here
"""
function build_sbin(
    i :: Int,
    j :: Int,
    gb :: BinomialSet{T, S}
) :: T where {T <: GBElement, S <: GBOrder}
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

function initial_gb(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    T :: DataType = Binomial,
) :: Vector{T}
    if T == Binomial
        gb = []
        #The current matrix A has n + m rows, 2n + m cols, so n = cols - rows
        num_gens = size(A, 2) - size(A, 1)
        for i in 1:num_gens
            gen = lattice_generator_binomial(i, A, b, C, u)
            if !isnothing(gen)
                push!(gb, gen)
            end
        end
    else
        gb = []
        num_gens = size(A, 2)
        for i in 1:num_gens
            gen = lattice_generator_graded(i, A, b, C, u)
            if !isnothing(gen)
                push!(gb, gen)
            end
        end
    end
    return gb
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
    minimization = structure == Binomial
    A, b, C, u = GBTools.normalize(
        A, b, C, u, apply_normalization=minimization
    )
    initial_basis = initial_gb(A, b, C, u, T = structure)
    gb = BinomialSet(initial_basis, MonomialOrder(C), minimization)
    i = 1
    iteration_count = 0
    spair_count = 0
    zero_reductions = 0
    #Main loop: generate all relevant S-binomials and reduce them
    while i <= length(gb)
        for j in 1:(i-1)
            iteration_count += 1
            if is_support_reducible(i, j, gb)
                continue
            end
            r = build_sbin(i, j, gb)
            if isfeasible(r, A, b, u)
                spair_count += 1
                reduced_to_zero = SupportTrees.reduce!(r, gb, reduction_tree(gb))
                if reduced_to_zero
                    zero_reductions += 1
                    continue
                end
                push!(gb, r)
            end
            if iteration_count % auto_reduce_freq == 0
                auto_reduce(gb)
            end
        end
        i += 1
    end
    @info "Buchberger: S-binomials reduced" iteration_count spair_count zero_reductions
    #Convert the basis to the same format 4ti2 uses
    minimal_basis!(gb)
    output_basis = BinomialSets.fourti2_form(gb)
    return output_basis
end

end
