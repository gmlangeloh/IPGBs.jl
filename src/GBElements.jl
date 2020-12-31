"""
This module defines all the binomial data structures used in my implementations
of Buchberger's algorithm and Signature-based algorithms.

TODO check whether this should be broken down into multiple modules

TODO there are multiple implementations of sparse subtraction of GradedBinomials
here. I should probably check whether I need them all
"""
module GBElements
export GBElement, regular_reducible,
    GradedBinomial, iszero, degree_reducible, supports, lattice_generator,
    isfeasible, fourti2_form, getfilter

using StaticArrays
using IPGBs.FastBitSets

"""
Gets the filter of a binomial, that is, the list of indices of variables
appearing in its leading term.

If fullfilter = true, include the indices of variables appearing in its
trailing term.
"""
function getfilter(
    binomial :: T;
    fullfilter :: Bool = false
) :: Vector{Int} where {T <: AbstractVector{Int}}
    filter = Int[]
    for i in 1:length(binomial)
        if binomial[i] > 0 || (fullfilter && binomial[i] != 0)
            push!(filter, i)
        end
    end
    return filter
end

"""
Checks whether `reducer` divides `g`, using the filter of `reducer` for
efficiency. When fullfilter = true, checks if g.head >= reducer.head and
g.tail >= reducer.tail coordinate-wise, while also checking a degree criterion.
"""
function reduces(
    g :: T,
    filter :: Vector{Int},
    reducer :: T,
    gb :: S;
    fullfilter :: Bool = true,
    negative :: Bool = false
) :: Bool where {T <: AbstractVector{Int}, S <: AbstractVector{T}}
    sign :: Int = negative ? -1 : 1
    if fullfilter
        for i in filter
            if sign * g[i] * reducer[i] < 0 #Different signs, doesn't divide
                return false
            elseif abs(g[i]) < abs(reducer[i]) #reducer doesn't divide
                return false
            end
        end
        #Checks the truncation criterion of Thomas and Weismantel
        if !degree_reducible(reducer, g, negative=negative)
            return false
        end
        #This is used in signature-based algorithms. Non-signature algorithms
        #simply return true here
        if !regular_reducible(reducer, g, gb)
            println("Singular reducer found")
            @show reducer reducer.signature
            @show g g.signature
            return false
        end
    else
        for i in filter
            if sign * g[i] < reducer[i]
                return false
            end
        end
    end
    return true
end

"""
Finds the maximum k such that k * reducer <= binomial coordinate-wise.
"""
function reduction_factor(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) :: Int where {T <: AbstractVector{Int}}
    i = 1
    while reducer[i] <= 0
        i += 1
    end
    factor = Int(floor(binomial[i] / reducer[i]))
    if (!negative && factor == 1) || (negative && factor == -1)
        return factor
    end
    i += 1
    while i <= length(reducer)
        if reducer[i] > 0
            newfactor = Int(floor(binomial[i] / reducer[i]))
            if (!negative && newfactor < factor) ||
                (negative && newfactor > factor)
                factor = newfactor
                if (!negative && factor == 1) || (negative && factor == -1)
                    return factor
                end
            end
        end
        i += 1
    end
    return factor
end

"""
Reduce `binomial` by `reducer` as many times as possible. Assumes `reducer`
reduces `binomial` at least once.
"""
function reduce!(
    binomial :: T,
    reducer :: T;
    negative :: Bool = false
) where {T <: AbstractVector{Int}}
    factor = reduction_factor(binomial, reducer, negative=negative)
    if !negative && factor == 1
        for i in 1:length(binomial)
            binomial[i] -= reducer[i]
        end
    elseif negative && factor == -1
        for i in 1:length(binomial)
            binomial[i] += reducer[i]
        end
    else
        for i in 1:length(binomial)
            binomial[i] -= factor * reducer[i]
        end
    end
end

"""
Abstract type used for GB computations. It is meant to generalize both binomials
and binomials with signature, in order to simplify the implementation of
reduction algorithms.
"""
abstract type GBElement <: AbstractVector{Int} end

"""
Returns true iff `reducer` is regular wrt `g`, that is, sig(reducer) < sig(g).
This is used when checking for regular reductions in a signature-based
algorithm. For a non-signature-based algorithm, and thus a generic GBElement,
this is always true.
"""
function regular_reducible(
    reducer :: GBElement,
    g :: GBElement,
    gb :: AbstractVector{T}
) :: Bool where {T <: GBElement}
    return true
end

"""
Represents a binomial as a vector of integers, tracking its leading and trailing
terms and degree for easy and efficient implementation of truncated GBs.
"""
mutable struct GradedBinomial <: GBElement
    element :: Vector{Int} #The element represented as a vector
    head :: Vector{Int} #The (ordered) indices of the positive coordinates of g
    tail :: Vector{Int} #The (ordered) indices of the negative coordinates of g
    cost :: Int #The cost of this element, for efficiency
    degree :: Vector{Int} #The degree of this element, for efficiency
    positive_degree :: Vector{Int} #Degree of the positive part of this element
    negative_degree :: Vector{Int} #Degree of the negative part of this element
end

function getfilter(
    binomial :: GradedBinomial;
    fullfilter :: Bool = false
) :: Vector{Int}
    if !fullfilter
        return binomial.head
    end
    #Return the filter ordered, merging both head and tail
    head_index = 1
    tail_index = 1
    head_unfinished = head_index <= length(binomial.head)
    tail_unfinished = tail_index <= length(binomial.tail)
    filter = Int[]
    while head_unfinished || tail_unfinished
        head_unfinished = head_index <= length(binomial.head)
        tail_unfinished = tail_index <= length(binomial.tail)
        if head_unfinished && tail_unfinished
            if binomial.head[head_index] < binomial.tail[tail_index]
                push!(filter, binomial.head[head_index])
                head_index += 1
                head_unfinished = head_index <= length(binomial.head)
            else
                push!(filter, binomial.tail[tail_index])
                tail_index += 1
                tail_unfinished = tail_index <= length(binomial.tail)
            end
        elseif head_unfinished #Tail already finished
            push!(filter, binomial.head[head_index])
            head_index += 1
            head_unfinished = head_index <= length(binomial.head)
        else #Head already finished
            push!(filter, binomial.tail[tail_index])
            tail_index += 1
            tail_unfinished = tail_index <= length(binomial.tail)
        end
    end
    return filter
end

#
# Implementation of the AbstractVector interface for GradedBinomials
#

function Base.size(
    g :: GradedBinomial
) :: Tuple
    return size(g.element)
end

function Base.getindex(
    g :: GradedBinomial,
    i :: Int
) :: Int
    return g.element[i]
end

"""
Sets `g` at index `i` to `v`. In theory, this shouldn't be used, as GradedBinomials
should not be changed.
"""
function Base.setindex!(
    g :: GradedBinomial,
    v :: Int,
    i :: Int
)
    g.element[i] = v
end

function Base.length(
    g :: GradedBinomial
) :: Int
    return length(g.element)
end

#
# Additional functions for a GradedBinomial
#

"""
Returns true iff g is smaller than h in the tiebreaking order (grevlex).
Assumes that g.cost == h.cost, that is, that there is a tie.

I don't really understand why this works, it doesn't look like grevlex to me.
But it gives precisely the same results as 4ti2, so I guess I'll keep it.
Commented below is an implementation which does not give the same results as
4ti2, but makes more sense to me.
"""
function lt_tiebreaker(
    g :: GradedBinomial,
    h :: GradedBinomial
) :: Bool
    @assert g.cost == h.cost
    gsmaller :: Int = 0 #-1 when g < h, 0 when g = h, 1 when g > h
    sum_g :: Int = 0
    sum_h :: Int = 0
    # Compute cumulative sums for g.element and h.element
    #the smallest one is the one with lowest cumulative sum at the farthest
    #point where they don't tie.
    for i in 1:length(g)
        sum_g += g[i]
        sum_h += h[i]
        if sum_g < sum_h
            gsmaller = -1
            break
        elseif sum_g > sum_h
            gsmaller = 1
            break
        end
    end
    return gsmaller == -1 ? true : false
    #for i in 1:length(g)
    #    sum_g += g[i]
    #    sum_h += h[i]
    #end
    #if sum_g < sum_h
    #    return true
    #elseif sum_g > sum_h
    #    return false
    #end
    #for i in 1:length(g)
    #    sum_g -= g[i]
    #    sum_h -= h[i]
    #    if sum_g < sum_h
    #        return true
    #    elseif sum_g > sum_h
    #        return false
    #    end
    #end
    #return false #If they are equal wrt grevlex at this point, g == h
end

"""
Decomposes g = v+ - v-, where both v+ and v- are non-negative arrays.
"""
function decomposition(
    g :: GradedBinomial
) :: Tuple{Vector{Int}, Vector{Int}}
    vp = zeros(Int, length(g.element))
    vm = zeros(Int, length(g.element))
    for i in g.head
        vp[i] = g.element[i]
    end
    for i in g.tail
        vm[i] = -g.element[i]
    end
    return vp, vm
end

"""
Given g = v+ - v-, returns the full form of g including slack variables.
This is given by:
x^(v+) s^((Av)-) r^(v-)   -   x^(v-) s^((Av)+) r^(v+)
"""
function fullform(
    g :: GradedBinomial
) :: Vector{Int}
    vp, vm = decomposition(g)
    Av = g.degree
    Avp = zeros(Int, length(Av))
    Avm = zeros(Int, length(Av))
    for i in 1:length(Av)
        if Av[i] > 0
            Avp[i] = Av[i]
        elseif Av[i] < 0
            Avm[i] = -Av[i]
        end
    end
    head = [ vp; Avm; vm ]
    tail = [ vm; Avp; vp ]
    return head - tail
end

"""
Checks whether g.degree+ <= h.degree coordinate-wise, where g.degree+ is
the set of positive subset of coordinates of g.degree.

If negative == true, checks whether g.degree- <= h.degree-
"""
function degree_reducible(
    g :: GradedBinomial,
    h :: GradedBinomial;
    negative :: Bool = false
) :: Bool
    m :: Int = length(g.degree)
    if !negative
        for i in 1:m
            if g.degree[i] > 0 && g.degree[i] > h.degree[i]
                return false
            end
        end
    else
        for i in 1:m
            #If both are negative and g is larger in abs value than h, then
            #g does not divide h
            if g.degree[i] < 0 && h.degree[i] < 0 && g.degree[i] < h.degree[i]
                return false
            end
        end
    end
    return true
end

function iszero(
    g :: GradedBinomial
) :: Bool
    return isempty(g.head)
end

function Base.show(
    io :: IO,
    g :: GradedBinomial
)
    print(io, g.element, " : c", g.cost, " : d", g.degree)
end

"""
Unary minus for GradedBinomials. It swaps head by tail and inverts cost.
"""
function Base.:-(
    g :: GradedBinomial
) :: GradedBinomial
    new_element = -g.element
    new_head = g.tail
    new_tail = g.head
    new_cost = -g.cost
    new_degree = -g.degree
    new_pos_degree = g.negative_degree
    new_neg_degree = g.positive_degree
    return GradedBinomial(new_element, new_head, new_tail, new_cost, new_degree,
                     new_pos_degree, new_neg_degree)
end

"""
In-place unary minus for GradedBinomials
"""
function opposite!(
    g :: GradedBinomial
) :: Nothing
    g.element .= .-(g.element)
		g.head, g.tail = g.tail, g.head
    g.cost = -g.cost
    g.degree .= .-(g.degree)
    g.positive_degree, g.negative_degree = g.negative_degree, g.positive_degree
		return
end

"""
Computes g = h - g changing g in-place. This is useful in the case g is
supposed to be reduced by h but g.cost < h.cost.
"""
function reduce_negative!(
    g :: GradedBinomial,
    h :: GradedBinomial
)
    g_h, g_t, h_h, h_t = (1, 1, 1, 1)
    while g_h <= length(g.head) || g_t <= length(g.tail) ||
        h_h <= length(h.head) || h_t <= length(h.tail)
        if g_h <= length(g.head)
            tgh = (g.head[g_h], 1)
        else
            tgh = (typemax(Int), 1)
        end
        if g_t <= length(g.tail)
            tgt = (g.tail[g_t], 2)
        else
            tgt = (typemax(Int), 2)
        end
        if h_h <= length(h.head)
            thh = (h.head[h_h], 3)
        else
            thh = (typemax(Int), 3)
        end
        if h_t <= length(h.tail)
            tht = (h.tail[h_t], 4)
        else
            tht = (typemax(Int), 4)
        end
        L = SVector(tgh, tgt, thh, tht)
        L = sort(L)
        min_idx = L[1][1]
        i = 1
        #First replace g[i] by -g[i], then add h[i]
        while i <= 4 && L[i][1] == min_idx
            if L[i][2] == 1 #g.head[g_h] == min_idx
                g.element[min_idx] = -g.element[min_idx]
            elseif L[i][2] == 2 #g.tail[g_t] == min_idx
                g.element[min_idx] = -g.element[min_idx]
            elseif L[i][2] == 3 #h.head[h_h] == min_idx
                g.element[min_idx] += h.element[min_idx]
                h_h += 1
            elseif L[i][2] == 4 #h.tail[h_t] == min_idx
                g.element[min_idx] += h.element[min_idx]
                h_t += 1
            end
            i += 1
        end
        #Update head / tail
        while g_h <= length(g.head) && g.head[g_h] < min_idx
            g_h += 1
        end
        while g_t <= length(g.tail) && g.tail[g_t] < min_idx
            g_t += 1
        end
        if g.element[min_idx] > 0 #Insert in head, remove from tail
            if g_h > length(g.head)
                push!(g.head, min_idx)
                g_h += 1
            elseif g.head[g_h] > min_idx
                insert!(g.head, g_h, min_idx)
                g_h += 1
            end
            if g_t <= length(g.tail) && g.tail[g_t] == min_idx
                deleteat!(g.tail, g_t)
                #g_t += 1
            end
        elseif g.element[min_idx] < 0 #Insert in tail, remove from head
            if g_t > length(g.tail)
                push!(g.tail, min_idx)
                g_t += 1
            elseif g.tail[g_t] > min_idx
                insert!(g.tail, g_t, min_idx)
                g_t += 1
            end
            if g_h <= length(g.head) && g.head[g_h] == min_idx
                deleteat!(g.head, g_h)
                #g_h += 1
            end
        else #g.element[min_idx] == 0 #Remove from head and tail
            if g_h <= length(g.head) && g.head[g_h] == min_idx
                deleteat!(g.head, g_h)
                #g_h += 1
            elseif g_t <= length(g.tail) && g.tail[g_t] == min_idx
                deleteat!(g.tail, g_t)
                #g_t += 1
            end
        end
    end
    g.cost = h.cost - g.cost
    for i in 1:length(g.degree)
        g.degree[i] = h.degree[i] - g.degree[i]
        #Yes, these are supposed to be swapped. The math works out.
        g.positive_degree[i] = h.negative_degree[i] - g.positive_degree[i]
        g.negative_degree[i] = h.positive_degree[i] - g.negative_degree[i]
    end
end

"""
Computes g -= h in place.
"""
function reduce!(
    g :: GradedBinomial,
    h :: GradedBinomial;
    negative :: Bool = false
)
    if negative || g.cost < h.cost || (g.cost == h.cost && lt_tiebreaker(g, h))
        reduce_negative!(g, h)
        return
    end
    g_h, g_t, h_h, h_t = (1, 1, 1, 1)
    while h_h <= length(h.head) || h_t <= length(h.tail)
        if h_h <= length(h.head)
            thh = (h.head[h_h], 3)
        else
            thh = (typemax(Int), 3)
        end
        if h_t <= length(h.tail)
            tht = (h.tail[h_t], 4)
        else
            tht = (typemax(Int), 4)
        end
        #TODO Don't even have to sort this with only 2 elements, could just use an if
        L = SVector(thh, tht)
        L = sort(L)
        min_idx = L[1][1]
        i = 1
        while i <= 2 && L[i][1] == min_idx
            if L[i][2] == 3
                g.element[min_idx] -= h.element[min_idx]
                h_h += 1
            elseif L[i][2] == 4
                g.element[min_idx] -= h.element[min_idx]
                h_t += 1
            end
            i += 1
        end
        #Update head / tail
        while g_h <= length(g.head) && g.head[g_h] < min_idx
            g_h += 1
        end
        while g_t <= length(g.tail) && g.tail[g_t] < min_idx
            g_t += 1
        end
        if g.element[min_idx] > 0 #Insert in head, remove from tail
            if g_h > length(g.head)
                push!(g.head, min_idx)
            elseif g.head[g_h] > min_idx
                insert!(g.head, g_h, min_idx)
            end
            if g_t <= length(g.tail) && g.tail[g_t] == min_idx
                deleteat!(g.tail, g_t)
            end
        elseif g.element[min_idx] < 0 #Insert in tail, remove from head
            if g_t > length(g.tail)
                push!(g.tail, min_idx)
            elseif g.tail[g_t] > min_idx
                insert!(g.tail, g_t, min_idx)
            end
            if g_h <= length(g.head) && g.head[g_h] == min_idx
                deleteat!(g.head, g_h)
            end
        else #g.element[min_idx] == 0 #Remove from head and tail
            if g_h <= length(g.head) && g.head[g_h] == min_idx
                deleteat!(g.head, g_h)
            elseif g_t <= length(g.tail) && g.tail[g_t] == min_idx
                deleteat!(g.tail, g_t)
            end
        end
    end
    g.cost -= h.cost
    for i in 1:length(g.degree)
        g.degree[i] -= h.degree[i]
        #Yes, these are supposed to be swapped. The math works out.
        g.positive_degree[i] -= h.negative_degree[i]
        g.negative_degree[i] -= h.positive_degree[i]
    end
end

"""
Subtracts two GradedBinomials, updating head/tail and cost.
This is more efficient than just subtracting directly.
"""
function Base.:-(
    g :: GradedBinomial,
    h :: GradedBinomial
) :: GradedBinomial
    new_element = zeros(Int, length(g.element))
    new_head = Int[]
    new_tail = Int[]
		# Henrique Becker Comment: Tested and only gives a slight slowdown,
		# new_head and new_tail probably are often much smaller.
		#sizehint!(new_head, length(g.element))
		#sizehint!(new_tail, length(g.element))
    g_h, g_t, h_h, h_t = (1, 1, 1, 1)
    while g_h <= length(g.head) || g_t <= length(g.tail) ||
        h_h <= length(h.head) || h_t <= length(h.tail)
        if g_h <= length(g.head)
            tgh = (g.head[g_h], 1)
        else
            tgh = (typemax(Int), 1)
        end
        if g_t <= length(g.tail)
            tgt = (g.tail[g_t], 2)
        else
            tgt = (typemax(Int), 2)
        end
        if h_h <= length(h.head)
            thh = (h.head[h_h], 3)
        else
            thh = (typemax(Int), 3)
        end
        if h_t <= length(h.tail)
            tht = (h.tail[h_t], 4)
        else
            tht = (typemax(Int), 4)
        end
        L = SVector(tgh, tgt, thh, tht)
        L = sort(L)
        min_idx = L[1][1]
        i = 1
        while i <= 4 && L[i][1] == min_idx
            if L[i][2] == 1
                new_element[min_idx] += g.element[min_idx]
                g_h += 1
            elseif L[i][2] == 2
                new_element[min_idx] += g.element[min_idx]
                g_t += 1
            elseif L[i][2] == 3
                new_element[min_idx] -= h.element[min_idx]
                h_h += 1
            elseif L[i][2] == 4
                new_element[min_idx] -= h.element[min_idx]
                h_t += 1
            end
            i += 1
        end
        if new_element[min_idx] > 0
            push!(new_head, min_idx)
        elseif new_element[min_idx] < 0
            push!(new_tail, min_idx)
        end
    end
    new_cost = g.cost - h.cost
    new_degree = g.degree - h.degree
    new_pos_degree = g.positive_degree - h.negative_degree
    new_neg_degree = g.negative_degree - h.positive_degree
    return GradedBinomial(new_element, new_head, new_tail, new_cost, new_degree,
                     new_pos_degree, new_neg_degree)
end

"""
Computes bitsets with the positive and negative supports of `g`.
"""
function supports(
    g :: GradedBinomial
) :: Tuple{FastBitSet, FastBitSet}
    pos_supp = Array(g.head)
    neg_supp = Array(g.tail)
    n = length(g)
    m = length(g.degree)
    for i in 1:m
        if g.degree[i] > 0
            push!(neg_supp, i + n)
        elseif g.degree[i] < 0
            push!(pos_supp, i + n)
        end
    end
    for i in g.head
        push!(neg_supp, i + n + m)
    end
    for i in g.tail
        push!(pos_supp, i + n + m)
    end
    bitset_length = 2n + m
    return makebitset(bitset_length, pos_supp), makebitset(bitset_length, neg_supp)
end

"""
Vector in Z^n with i-th coordinate 1 and remaining coordinates 0.
"""
function lattice_generator(
    i :: Int,
    A :: Array{Int, 2},
    c :: Array{Int}
) :: GradedBinomial
    n = size(A, 2)
    v = zeros(Int, n)
    v[i] = 1
    if ndims(c) == 1
        cost = c[i]
    else
        cost = c[1, i]
    end
    degree = A * v
    pos_degree = A * v
    neg_degree = zeros(Int, length(pos_degree))
    return GradedBinomial(v, Int[i], Int[], cost, degree, pos_degree, neg_degree)
end

"""
Returns true iff v1 <= v2 coordinate-wise.
"""
function le_coordinatewise(
    v1 :: Vector{Int},
    v2 :: Vector{Int}
) :: Bool
    return all(i -> v1[i] <= v2[i], keys(v1))
end

"""
Checks whether element is smaller coordinate-wise than u only considering the
given indices.
"""
function sparse_le(
    indices :: Vector{Int},
    element :: Vector{Int},
    u :: Vector{Int}
) :: Bool
    return all(i -> abs(element[i]) <= u[i], indices)
end

"""
Returns true iff v should be considered for reduction in a truncated GB
algorithm.
"""
function isfeasible(
    v :: GradedBinomial,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Bool
    head = v.positive_degree
    tail = v.negative_degree
    return le_coordinatewise(head, b) && le_coordinatewise(tail, b) &&
        sparse_le(v.head, v.element, u) && sparse_le(v.tail, v.element, u)
end

"""
Transforms a vector of GradedBinomials into the form of 4ti2 input/output.
This includes putting GradedBinomials into fullform, including all of their
variables which are implicit in the GradedBinomial structure.

If `is_maximization` is true, assumes the GradedBinomials are in maximization
form. As 4ti2 assumes minimization, signs have to be swapped.
"""
function fourti2_form(
    binomial_set :: Vector{GradedBinomial};
    is_maximization :: Bool = true
) :: Array{Int, 2}
    fourti2_elems = []
    for binomial in binomial_set
        f = GBElements.fullform(binomial)
        if is_maximization
            f = -f
        end
        push!(fourti2_elems, f)
    end
    fourti2_set = foldl(hcat, fourti2_elems)
    return fourti2_set'
end

end
