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
using IPGBs.Statistics
using IPGBs.SupportTrees

using IPGBs.GBAlgorithms

#4ti2 uses 2500 here. I think auto-reducing makes things slower, though...
const AUTO_REDUCE_FREQ = 100000

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
    removed_by_autoreduction :: Int

    #These fields may be specific to the Buchberger algorithm
    eliminated_by_gcd :: Int

    BuchbergerStats() = new(0, 0, 0, 0, 0, 0, 0, 0)
end

mutable struct BuchbergerAlgorithm{T <: GBElement} <: GBAlgorithm
    basis :: BinomialSet{T, MonomialOrder}
    state :: BuchbergerState
    should_truncate :: Bool
    stats :: BuchbergerStats
    preallocated :: Vector{Int}
    num_iterations :: Int #Total number of binomials added to the basis at some point

    function BuchbergerAlgorithm(
        T :: Type,
        C :: Array{Float64, 2},
        A :: Array{Int, 2},
        b :: Vector{Int},
        should_truncate :: Bool,
        minimization :: Bool
    )
        order = MonomialOrder(C, A, b)
        state = BuchbergerState(0)
        stats = BuchbergerStats()
        new{T}(
            BinomialSet{T, MonomialOrder}(T[], order, minimization), state,
            should_truncate, stats, Int[], 0
        )
    end
end

GBAlgorithms.use_implicit_representation(:: BuchbergerAlgorithm{T}) where {T} = is_implicit(T)

"""
Produces the next BinomialPair to be processed by the Buchberger algorithm,
if one exists, or nothing otherwise.

Auto-reduces the partial GB periodically (for consistency with 4ti2).
"""
function GBAlgorithms.next_pair!(
    algorithm :: BuchbergerAlgorithm{T}
) :: Union{BinomialPair, Nothing} where {T <: GBElement}
    previous_i = algorithm.state.i
    s = next_state!(algorithm.state)
    if !isnothing(s)
        i, j = s
        if i != previous_i #Started processing S-pairs generated by new binomial
            algorithm.num_iterations += 1
            #Auto-reduce the basis periodically
            if algorithm.num_iterations % AUTO_REDUCE_FREQ == 0 &&
                algorithm.num_iterations != 0
                removed, before_idx = auto_reduce_once!(current_basis(algorithm), i)
                algorithm.stats.removed_by_autoreduction += removed
                algorithm.state.n -= removed
                algorithm.state.i -= before_idx
                i -= before_idx
            end
        end
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
    push!(current_basis(algorithm), copy(g))
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
    C :: Array{Float64, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    change_ordering!(current_basis(algorithm), C, A, b)
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    num_vars = size(A, 2)
    algorithm.preallocated = Vector{Int}(undef, num_vars)
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u, order(algorithm.basis),
                              check_truncation=truncate_basis(algorithm))
        if !isnothing(e)
            update!(algorithm, e, nothing)
        end
    end
end

end
