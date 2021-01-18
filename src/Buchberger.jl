"""
Implements a combinatorial Buchberger algorithm for Integer Programming
where all data is non-negative. Based on Thomas and Weismantel (1997).
"""
module Buchberger
export BuchbergerAlgorithm

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.SupportTrees

using IPGBs.GBAlgorithms

"""
The state of Buchberger S-binomial generation.
"""
mutable struct BuchbergerState
    i :: Int
    j :: Int
    n :: Int
    BuchbergerState(n) = new(1, 1, n)
end

increment_size!(state :: BuchbergerState) = state.n += 1

"""
Updates the state to the indices of the next generated S-pair. Returns the indices
of the new state, if it corresponds to a new S-pair, or `nothing`, if all S-pairs
were already generated.
"""
function next_state!(
    state :: BuchbergerState
) :: Union{Tuple{Int, Int}, Nothing}
    if state.j < state.i - 1
        state.j += 1
        return (state.i, state.j)
    elseif state.i < state.n
        state.i += 1
        state.j = 1
        return (state.i, state.j)
    end
    return nothing
end

struct BuchbergerAlgorithm{T <: GBElement} <: GBAlgorithm
    basis :: BinomialSet{T, MonomialOrder}
    state :: BuchbergerState
    stats :: GBStats

    function BuchbergerAlgorithm(T :: Type, C :: Array{Int, 2})
        order = MonomialOrder(C)
        state = BuchbergerState(0)
        stats = GBStats()
        stats.stats["eliminated_by_gcd"] = 0
        new{T}(BinomialSet{T, MonomialOrder}(T[], order), state, stats)
    end
end

stats(algorithm :: BuchbergerAlgorithm) = algorithm.stats.stats
current_basis(algorithm :: BuchbergerAlgorithm) = algorithm.basis

function next_pair!(
    algorithm :: BuchbergerAlgorithm{T}
) :: Union{BinomialPair, Nothing} where {T <: GBElement}
    s = next_state!(algorithm.state)
    if !isnothing(s)
        i, j = s
        stats(algorithm)["queued_pairs"] += 1
        return BinomialPair(i, j)
    end
    return nothing
end

function update!(
    algorithm :: BuchbergerAlgorithm{T},
    g :: T
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    increment_size!(algorithm.state)
    stats(algorithm)["max_basis_size"] = max(stats(algorithm)["max_basis_size"],
                                             length(current_basis(algorithm)))
end

"""
Applies the GCD criterion to determine whether or not to eliminate the given
S-pair.
"""
function late_pair_elimination(
    algorithm :: BuchbergerAlgorithm{T},
    pair :: CriticalPair
) :: Bool where {T <: GBElement}
    if is_support_reducible(first(pair), second(pair), current_basis(algorithm))
        stats(algorithm)["eliminated_by_gcd"] += 1
        return true
    end
    return false
end

function process_zero_reduction!(
    algorithm :: BuchbergerAlgorithm{T},
    :: T
) where {T <: GBElement}
    stats(algorithm)["zero_reductions"] += 1
end

function initialize!(
    algorithm :: BuchbergerAlgorithm{T},
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    num_vars = size(A, 2)
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u)
        if !isnothing(e)
            update!(algorithm, e)
        end
    end
end

#TODO refactor from this point onward

#function initial_gb(
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    C :: Array{Int, 2},
#    u :: Vector{Int};
#    T :: DataType = Binomial,
#) :: Vector{T}
#    if T == Binomial
#        gb = []
#        #The current matrix A has n + m rows, 2n + m cols, so n = cols - rows
#        num_gens = size(A, 2) - size(A, 1)
#        for i in 1:num_gens
#            gen = lattice_generator_binomial(i, A, b, C, u)
#            if !isnothing(gen)
#                push!(gb, gen)
#            end
#        end
#    else
#        gb = []
#        num_gens = size(A, 2)
#        for i in 1:num_gens
#            gen = lattice_generator_graded(i, A, b, C, u)
#            if !isnothing(gen)
#                push!(gb, gen)
#            end
#        end
#    end
#    return gb
#end
#
#"""
#Computes a test set / Gröbner Basis for the IP:
#max C^T * x
#s.t. A * x <= b
#x <= u
#
#Structure refers to the data structure used to represent binomials internally.
#It can be either `Binomial` or `GradedBinomial`.
#"""
#function buchberger(
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    C :: Array{Int, 2},
#    u :: Vector{Int};
#    structure :: DataType = Binomial,
#    auto_reduce_freq :: Int = 2500
#) :: Vector{Vector{Int}}
#    @assert structure == Binomial || structure == GradedBinomial
#    minimization = structure == Binomial
#    A, b, C, u = GBTools.normalize(
#        A, b, C, u, apply_normalization=minimization
#    )
#    initial_basis = initial_gb(A, b, C, u, T = structure)
#    gb = BinomialSet(initial_basis, MonomialOrder(C), minimization)
#    i = 1
#    iteration_count = 0
#    spair_count = 0
#    zero_reductions = 0
#    #Main loop: generate all relevant S-binomials and reduce them
#    while i <= length(gb)
#        for j in 1:(i-1)
#            iteration_count += 1
#            if is_support_reducible(i, j, gb)
#                continue
#            end
#            r = build_sbin(i, j, gb)
#            if isfeasible(r, A, b, u)
#                spair_count += 1
#                reduced_to_zero = BinomialSets.reduce!(r, gb)
#                if reduced_to_zero
#                    zero_reductions += 1
#                    continue
#                end
#                push!(gb, r)
#            end
#            if iteration_count % auto_reduce_freq == 0
#                auto_reduce!(gb)
#            end
#        end
#        i += 1
#    end
#    @info "Buchberger: S-binomials reduced" iteration_count spair_count zero_reductions
#    #Convert the basis to the same format 4ti2 uses
#    minimal_basis!(gb)
#    output_basis = BinomialSets.fourti2_form(gb)
#    return output_basis
#end

end
