module Pairs

using DataStructures

abstract type BuchbergerState end

increment_size(:: BuchbergerState) = error("Not implemented.")
next_state!(:: BuchbergerState) = error("Not implemented.")

mutable struct PriorityState <: BuchbergerState
    heap :: BinaryMinHeap{Tuple{Int, Int}}
    n :: Int

    function PriorityState(n)
        #TODO: use Base.by to define how tuples should be compared
        heap = BinaryMinHeap{Tuple{Int, Int}}()
        for i in 1:n
            for j in 1:i-1
                push!(heap, (i, j))
            end
        end
        new(heap, n)
    end
end

mutable struct BatchState <: BuchbergerState
    processed :: Vector{Bool}
    queue :: BinaryMinHeap{Int}
    j :: Int
    n :: Int
    BatchState(n) = new(Bool[], n)
end

mutable struct FIFOState <: BuchbergerState
    i :: Int
    j :: Int
    n :: Int
    FIFOState(n) = new(1, 1, n)
end

increment_size!(state :: FIFOState) = state.n += 1

"""
Updates the state to the indices of the next generated S-pair. Returns the indices
of the new state, if it corresponds to a new S-pair, or `nothing`, if all S-pairs
were already generated.
"""
function next_state!(
    state :: FIFOState
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

end
