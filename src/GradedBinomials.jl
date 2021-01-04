"""
A structure for efficient computation of truncated GBs of IPs where the data
is positive, as described in Section 3 of Thomas and Weismantel (1997).

TODO there are multiple implementations of sparse subtraction of GradedBinomials
here. I should probably check whether I need them all
"""
module GradedBinomials
export GradedBinomial, lattice_generator, isfeasible, fourti2_form, supports

using StaticArrays
using IPGBs.FastBitSets
using IPGBS.GBElements

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
function IPGBs.GBElements.degree_reducible(
    g :: GradedBinomial,
    h :: GradedBinomial,
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
function IPGBs.GBElements.opposite!(
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
function IPGBs.GBElements.reduce!(
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
function IPGBs.GBElements.supports(
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
function IPGBs.GBElements.lattice_generator(
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

function degrees(
    v :: GradedBinomial,
    A :: Array{Int, 2}
) :: Tuple{Vector{Int}, Vector{Int}}
    return v.positive_degree, v.negative_degree
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
