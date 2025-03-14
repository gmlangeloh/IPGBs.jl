module Pairs

using DataStructures

using IPGBs.BinomialSets

export BuchbergerState, initialize_pairs

abstract type BuchbergerState end
abstract type PriorityState <: BuchbergerState end

update!(:: BuchbergerState) = error("Not implemented.")
next_pair!(:: BuchbergerState) = error("Not implemented.")
auto_reduce_now(:: BuchbergerState) = false
remove_auto_reduced!(:: BuchbergerState, :: Int) = error("Not implemented.")

function initialize_pairs(
    n :: Int,
    strategy :: Symbol,
    gb :: BinomialSet{T},
    x :: Union{T, Nothing}
) where {T <: AbstractVector{Int}}
    if strategy == :Pair
        return PairState(n, gb, x)
    elseif strategy == :Batch
        return BatchState(n, gb, x)
    elseif strategy == :FIFO
        return FIFOState(n)
    else
        error("Unknown strategy: $strategy")
    end
end

heap(:: PriorityState) = error("Not implemented.")
size(:: PriorityState) = error("Not implemented.")
increment_size!(:: PriorityState) = error("Not implemented.")

function update!(state :: PriorityState)
    increment_size!(state)
    n = size(state)
    for j in 1:(n - 1)
        push!(heap(state), (n, j))
    end
end

function next_pair!(
    state :: PriorityState
) :: Tuple{Int, Int}
    if isempty(heap(state))
        return (-1, -1)
    end
    return pop!(heap(state))
end

function by_solution(
    i :: Int,
    gb :: BinomialSet{T},
    x :: T
) where {T <: AbstractVector{Int}}
    g = gb[i]
    y = x + g
    if all(yi >= 0 for yi in y)
        return (0, i)
    end
    return (1, i)
end

function by_pair_solution(
    i :: Int,
    j :: Int,
    gb :: BinomialSet{T},
    x :: T
) where {T <: AbstractVector{Int}}
    a1, _ = by_solution(i, gb, x)
    a2, _ = by_solution(j, gb, x)
    if a1 == 1 && a2 == 1
        return (1, i, j)
    end
    return (0, i, j)
end

mutable struct PairState <: PriorityState
    heap :: BinaryHeap{Tuple{Int, Int}}
    n :: Int

    function PairState(n, gb, x)
        heap = BinaryMinHeap{Tuple{Int, Int}}()
        if !isnothing(x)
            f((i, j)) = by_pair_solution(i, j, gb, x)
            heap = BinaryHeap{Tuple{Int, Int}}(Base.By(f))
        end
        for i in 1:n
            for j in 1:i-1
                push!(heap, (i, j))
            end
        end
        new(heap, n)
    end
end

heap(state :: PairState) = state.heap
size(state :: PairState) = state.n
increment_size!(state :: PairState) = state.n += 1

mutable struct BatchState <: PriorityState
    heap :: BinaryHeap{Int}
    current_i :: Int
    j :: Int
    n :: Int

    function BatchState(
        n :: Int,
        gb :: BinomialSet{T},
        x :: Union{T, Nothing}
    ) where {T <: AbstractVector{Int}}
        heap = BinaryMinHeap{Int}(1:n)
        if !isnothing(x)
            f(i) = by_solution(i, gb, x)
            heap = BinaryHeap{Int}(Base.By(f), 1:n)
        end
        new(heap, 0, 1, n)
    end
end

heap(state :: BatchState) = state.heap
size(state :: BatchState) = state.n
increment_size!(state :: BatchState) = state.n += 1

function update!(state :: BatchState)
    increment_size!(state)
    push!(heap(state), size(state))
end

function next_pair!(state :: BatchState) :: Tuple{Int, Int}
    if state.j < state.current_i - 1
        state.j += 1
        return (state.current_i, state.j)
    elseif isempty(heap(state))
        return nothing
    else # state.j == state.current_i, get a new one
        new_i = pop!(heap(state))
        if new_i == 1 # There is no pair with i == 1, as j < i
            new_i = pop!(heap(state))
        end
        state.current_i = new_i
        state.j = 1
        return (state.current_i, state.j)
    end
    return (-1, -1)
end

mutable struct FIFOState <: BuchbergerState
    i :: Int
    j :: Int
    n :: Int
    FIFOState(n) = new(1, 1, n)
end

increment_size!(state :: FIFOState) = state.n += 1
update!(state :: FIFOState) = increment_size!(state)

"""
Updates the state to the indices of the next generated S-pair. Returns the indices
of the new state, if it corresponds to a new S-pair, or `nothing`, if all S-pairs
were already generated.
"""
function next_pair!(
    state :: FIFOState
) :: Tuple{Int, Int}
    if state.j < state.i - 1
        state.j += 1
        return (state.i, state.j)
    elseif state.i < state.n
        state.i += 1
        state.j = 1
        return (state.i, state.j)
    end
    return (-1, -1)
end

auto_reduce_now(state :: FIFOState) = state.j == state.i - 1

function remove_auto_reduced!(
    state :: FIFOState,
    removed :: Int,
    before_i :: Int
)
    state.n -= removed
    state.i -= before_i
end

end
