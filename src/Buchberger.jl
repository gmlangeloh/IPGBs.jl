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
using IPGBs.Orders
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

mutable struct BuchbergerStats <: GBStats
    #These fields will be present in any GBStats type
    zero_reductions :: Int
    max_basis_size :: Int
    queued_pairs :: Int
    built_pairs :: Int
    reduced_pairs :: Int
    eliminated_by_truncation :: Int

    #These fields may be specific to the Buchberger algorithm
    eliminated_by_gcd :: Int

    BuchbergerStats() = new(0, 0, 0, 0, 0, 0, 0)
end

struct BuchbergerAlgorithm{T <: GBElement} <: GBAlgorithm
    basis :: BinomialSet{T, MonomialOrder}
    state :: BuchbergerState
    should_truncate :: Bool
    stats :: BuchbergerStats

    function BuchbergerAlgorithm(
        T :: Type,
        C :: Array{Int, 2},
        should_truncate :: Bool,
        minimization :: Bool
    )
        order = MonomialOrder(C)
        state = BuchbergerState(0)
        stats = BuchbergerStats()
        new{T}(
            BinomialSet{T, MonomialOrder}(T[], order, minimization), state,
            should_truncate, stats
        )
    end
end

GBAlgorithms.use_implicit_representation(:: BuchbergerAlgorithm{T}) where {T} = is_implicit(T)

function GBAlgorithms.next_pair!(
    algorithm :: BuchbergerAlgorithm{T}
) :: Union{BinomialPair, Nothing} where {T <: GBElement}
    s = next_state!(algorithm.state)
    if !isnothing(s)
        i, j = s
        increment(algorithm, :queued_pairs)
        return BinomialPair(i, j)
    end
    return nothing
end

function GBAlgorithms.update!(
    algorithm :: BuchbergerAlgorithm{T},
    g :: T,
    :: Union{CriticalPair, Nothing} = nothing
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    increment_size!(algorithm.state)
    algorithm.stats.max_basis_size = max(algorithm.stats.max_basis_size,
                                         length(current_basis(algorithm)))
end

"""
Applies the GCD criterion to determine whether or not to eliminate the given
S-pair.
"""
function GBAlgorithms.late_pair_elimination(
    algorithm :: BuchbergerAlgorithm{T},
    pair :: CriticalPair
) :: Bool where {T <: GBElement}
    if is_support_reducible(GBElements.first(pair), GBElements.second(pair),
                            current_basis(algorithm))
        increment(algorithm, :eliminated_by_gcd)
        return true
    end
    return false
end

function GBAlgorithms.process_zero_reduction!(
    algorithm :: BuchbergerAlgorithm{T},
    :: T,
    :: CriticalPair
) where {T <: GBElement}
    increment(algorithm, :zero_reductions)
end

function GBAlgorithms.initialize!(
    algorithm :: BuchbergerAlgorithm{T},
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    change_ordering!(current_basis(algorithm), C) #TODO maybe this should be done in GBAlgorithm?
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    num_vars = size(A, 2)
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u, order(algorithm.basis),
                              check_truncation=truncate_basis(algorithm))
        if !isnothing(e)
            update!(algorithm, e, nothing)
        end
    end
end

end
