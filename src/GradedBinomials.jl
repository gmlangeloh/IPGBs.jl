"""
A structure for efficient computation of truncated GBs of IPs where the data
is positive, as described in Section 3 of Thomas and Weismantel (1997).
This type of implicit representation of GB elements is also used in Urbaniak,
Weismantel and Ziegler (1997).

TODO there are multiple implementations of sparse subtraction of GradedBinomials
here. I should probably check whether I need them all

TODO Currently, using GradedBinomials gives different (but apparently correct)
GBs than using Binomials. There's probably something implicit about the monomial
order they're computed wrt that I'm missing, due to one being maximization and the
other being minimization.

Regardless, I'm unconvinced this will ever be more efficient than the usual Binomials.
The main reason is that the implicit variables are updated and checked in the
reduction process regardless, only this happens in a different way in the SupportTree.
A slightly different implementation, where the implicit variables aren't stored
anywhere and are computed as needed could be useful if memory became an issue,
but that is another thing entirely.
"""
module GradedBinomials
export GradedBinomial, lattice_generator_graded, fourti2_form

using StaticArrays
using IPGBs.FastBitSets
using IPGBs.GBElements
using IPGBs.Orders

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

function Base.copy(
    g :: GradedBinomial
) :: GradedBinomial
    #TODO Check if I really only need to copy g.element.
    #It is quite possible I need to copy at least the degrees as well...
    return GradedBinomial(
        copy(g.element),
        g.head,
        g.tail,
        g.cost,
        g.degree,
        g.positive_degree,
        g.negative_degree
    )
end

"""
TODO do this in a memory-efficient way by reusing `result`
"""
function GBElements.minus(
    result :: Vector{Int},
    g :: GradedBinomial,
    h :: GradedBinomial
) :: GradedBinomial
    return g - h
end

GBElements.is_implicit(:: Type{<: GradedBinomial}) = true

"""
Computes the filter of `binomial`. In case fullfilter == false, this is just
the indices in the support of the leading term of `binomial`.

Otherwise (fullfilter == true) this is the vector of ordered indices of the
support of `binomial`, so it also includes indices in the trailing term.
"""
function filter(
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

GBElements.cost(g :: GradedBinomial) = g.cost

"""
Given g = v+ - v-, returns the full form of g including slack variables.
This is given by:
x^(v+) s^((Av)-) r^(v-)   -   x^(v-) s^((Av)+) r^(v+)
"""
function GBElements.fullform(
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
function GBElements.degree_reducible(
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

function GBElements.leading_term(
    g :: GradedBinomial
) :: Vector{Int}
    lt = zeros(Int, length(g))
    for i in g.head
        lt[i] = g[i]
    end
    return lt
end

GBElements.head(g :: GradedBinomial) = g.head
GBElements.is_zero(g :: GradedBinomial) = isempty(head(g))

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
function GBElements.opposite!(
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
) :: Bool
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
    return GBElements.is_zero(g)
end

"""
Computes g -= h in place.

Returns true iff g reduced to zero.
"""
function GBElements.reduce!(
    g :: GradedBinomial,
    h :: GradedBinomial,
    order :: GBOrder;
    negative :: Bool = false
) :: Bool
    if negative || Base.lt(order, g, h)
        return reduce_negative!(g, h)
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
                g[min_idx] -= h.element[min_idx]
                if g[min_idx] != 0
                end
                h_h += 1
            elseif L[i][2] == 4
                g[min_idx] -= h.element[min_idx]
                if g[min_idx] != 0
                end
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
        if g[min_idx] > 0 #Insert in head, remove from tail
            if g_h > length(g.head)
                push!(g.head, min_idx)
            elseif g.head[g_h] > min_idx
                insert!(g.head, g_h, min_idx)
            end
            if g_t <= length(g.tail) && g.tail[g_t] == min_idx
                deleteat!(g.tail, g_t)
            end
        elseif g[min_idx] < 0 #Insert in tail, remove from head
            if g_t > length(g.tail)
                push!(g.tail, min_idx)
            elseif g.tail[g_t] > min_idx
                insert!(g.tail, g_t, min_idx)
            end
            if g_h <= length(g.head) && g.head[g_h] == min_idx
                deleteat!(g.head, g_h)
            end
        else #g[min_idx] == 0 #Remove from head and tail
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
    return GBElements.is_zero(g)
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
function GBElements.supports(
    g :: GradedBinomial
) :: Tuple{FastBitSet, FastBitSet}
    pos_supp = Vector(g.head)
    neg_supp = Vector(g.tail)
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
    return FastBitSet(bitset_length, pos_supp), FastBitSet(bitset_length, neg_supp)
end

"""
Vector in Z^n with i-th coordinate 1 and remaining coordinates 0.
"""
function lattice_generator_graded(
    i :: Int,
    A :: Array{Int, 2},
    b :: Vector{Int},
    c :: Array{Float64},
    u :: Vector{Int},
    :: Union{GBOrder, Nothing} = nothing;
    check_truncation :: Bool = true
) :: Union{GradedBinomial, Nothing}
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
    generator = GradedBinomial(
        v, Int[i], Int[], cost, degree, pos_degree, neg_degree
    )
    if !check_truncation || isfeasible(generator, A, b, u)
        return generator
    end
    return nothing
end

function GBElements.degrees(
    v :: GradedBinomial,
    :: Array{Int, 2}
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
        f = fullform(binomial)
        if is_maximization
            f = -f
        end
        push!(fourti2_elems, f)
    end
    fourti2_set = foldl(hcat, fourti2_elems)
    return fourti2_set'
end

end
