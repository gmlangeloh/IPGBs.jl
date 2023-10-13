"""
Implements a combinatorial Buchberger algorithm for Integer Programming
where all data is non-negative. Based on Thomas and Weismantel (1997).
"""
module Buchberger
export BuchbergerAlgorithm

using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances
using MIPMatrixTools.SolverTools

using IPGBs
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.Orders
using IPGBs.Statistics
using IPGBs.SupportTrees

using IPGBs.GBAlgorithms

using JuMP

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
    truncation_type :: Symbol
    stats :: BuchbergerStats
    preallocated :: Vector{Int}
    num_iterations :: Int #Total number of binomials added to the basis at some point
    model :: JuMP.Model
    model_vars :: Vector{VariableRef}
    model_constraints :: Vector{ConstraintRef}
    instance :: IPInstance
    truncated_gens :: Vector{T}
    truncation_weight :: Vector{Float64}
    max_truncation_weight :: Float64

    function BuchbergerAlgorithm(
        markov :: Vector{Vector{Int}},
        instance :: IPInstance;
        T :: Type = Binomial,
        minimization :: Bool = true,
        truncation_type :: Symbol = :Model,
        trunc_var_type :: DataType = Real
    )
        @assert !isempty(markov)
        #Build order and generating set
        order = MonomialOrder(
            instance.C, instance.A, instance.b, 
            unbounded_variables(instance), minimization
        )
        initialize_binomials(instance, order)
        generating_set = [to_gbelement(m, order, T) for m in markov]
        preallocated = Vector{Int}(undef, length(generating_set[1]))
        should_truncate = truncation_type != :None
        #Initialize a feasibility model in case we want to use model truncation
        model, vars, constrs = SolverTools.feasibility_model(
            instance.A, instance.b, instance.u, nonnegative_vars(instance),
            trunc_var_type
        )
        #Truncate elements in the generating set
        nontruncated_gens = T[]
        truncated_gens = T[]
        if should_truncate
            @debug "Truncating generating set with algorithm: $(truncation_type)"
        end
        for gen in generating_set
            if is_zero(gen)
               nontruncated_gens = [gen] 
               truncated_gens = T[]
               break
            end
            truncated = GBElements.truncate(
                gen, instance.A, instance.b, instance.u, model, constrs,
                should_truncate, truncation_type
            )
            if truncated
                push!(truncated_gens, gen)
            else
                push!(nontruncated_gens, gen)
            end
        end
        binomial_gen_set = BinomialSet{T, MonomialOrder}(
            nontruncated_gens, order, minimization
        )
        auto_reduce_once!(binomial_gen_set)
        weight, max_weight = truncation_weight(instance)
        #Initialize the state of the algorithm (no pairs processed yet)
        state = BuchbergerState(length(binomial_gen_set))
        stats = BuchbergerStats()
        new{T}(
            binomial_gen_set, state, should_truncate, truncation_type, stats, 
            preallocated, 0, model, vars, constrs, instance, truncated_gens, 
            weight, max_weight
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
    #If we're starting to produce S-binomials from a new basis element
    #then check whether we should auto-reduce the basis for efficiency
    if algorithm.state.j == algorithm.state.i - 1
        algorithm.num_iterations += 1
        if IPGBs.should_auto_reduce(algorithm.num_iterations, length(current_basis(algorithm)))
            removed, before_idx = auto_reduce_once!(current_basis(algorithm), current_index=algorithm.state.i+1)
            algorithm.stats.removed_by_autoreduction += removed
            algorithm.state.n -= removed
            algorithm.state.i -= before_idx
        end
    end
    #Produce the next S-pair itself
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

"""
    GBAlgorithms.quick_truncation(
    algorithm :: BuchbergerAlgorithm{T},
    g :: T
) where {T <: GBElement}

Return true iff `g` is truncated by the weight criterion, that is, its
weight is higher than the maximum weight for fibers of the given instance.
"""
function GBAlgorithms.quick_truncation(
    algorithm :: BuchbergerAlgorithm{T},
    g :: T
) where {T <: GBElement}
    return weight(g, algorithm.truncation_weight) > algorithm.max_truncation_weight
end

end
