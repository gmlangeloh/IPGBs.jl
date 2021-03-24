"""
Implementation of the basic signature-polynomial pairs structures.

TODO this whole module would probably be cleaner if it was separated in
multiple modules. One defining signatures, another defining module monomial
orders, another defining SigPolys and so on.
"""
module SignaturePolynomials
export Signature, SigPoly, SigBasis, ModuleMonomialOrdering,
    pot_order, ltpot_order, top_order, is_zero, koszul,
    signature, SignaturePair, regular_spair, SigLead

using IPGBs.BinomialSets
using IPGBs.FastBitSets
using IPGBs.FastComparator
using IPGBs.GBElements
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.Orders
using IPGBs.SupportTrees

#
# Signature implementation
# They are a type of AbstractVector in order to use SupportTrees for divisor
# queries.
#

"""
A signature, or module monomial, is represented by a monomial with an index,
where the monomial is the coefficient of the module basis vector of the given
index.

TODO I'm thinking that `monomial` could easily be a sparse vector in cases with
many variables. I expect many entries to be 0 in this case, but it is worth it
to experiment with that.
"""
struct Signature <: AbstractVector{Int}
    index :: Int
    monomial :: Vector{Int}
end

function Base.isequal(
    s1 :: Signature,
    s2 :: Signature
) :: Bool
    return s1.index == s2.index && s1.monomial == s2.monomial
end

function Base.show(
    io :: IO,
    s :: Signature
)
    print(io, s.monomial, "e", s.index)
end

"""
Returns true iff `s1` divides `s2`.
"""
function divides(
    s1 :: Signature,
    s2 :: Signature
) :: Bool
    if s1.index != s2.index
        return false
    end
    n = length(s1.monomial)
    @assert n == length(s2.monomial)
    return all(i -> s1.monomial[i] <= s2.monomial[i], 1:n)
end

function Base.:*(
    monomial :: Vector{Int},
    s :: Signature
) :: Signature
    new_monomial = Vector{Int}(undef, length(monomial))
    @assert length(new_monomial) == length(s.monomial)
    for i in 1:length(monomial)
        new_monomial[i] = monomial[i] + s.monomial[i]
    end
    return Signature(s.index, new_monomial)
end

function Base.size(
    s :: Signature
) :: Tuple
    return size(s.monomial)
end

function Base.getindex(
    s :: Signature,
    i :: Int
) :: Int
    return s.monomial[i]
end

function Base.setindex!(
    s :: Signature,
    v :: Int,
    i :: Int
)
    s.monomial[i] = v
end

function Base.length(
    s :: Signature
) :: Int
    return length(s.monomial)
end

#
# Sig-lead ratios: represented as Signatures, but the monomial field may have
# some negative coordinates. SigLead comparison is implemented in the exact same
# way than signature comparison.
#

const SigLead = Signature

"""
Computes the SigLead ratio of `binomial` with `signature`. This is defined as
a signature with the same index as `signature`, but with monomial defined by
`signature.monomial - leading_term(binomial)`
"""
function siglead(
    binomial :: T,
    signature :: Signature
) :: SigLead where {T <: GBElement}
    ratio = signature.monomial - leading_term(binomial)
    return SigLead(signature.index, ratio)
end

"""
Represents a polynomial with its signature. Used in signature-based algorithms.
Stores its sig-lead ratio to compare signatures more efficiently when building
S-pairs, Koszul syzygies, reductions and elimination criteria.
"""
struct SigPoly{T <: GBElement} <: GBElement
    polynomial :: T
    signature :: Signature
    siglead :: SigLead
    SigPoly{T}(p :: T, s :: Signature) where {T <: GBElement} =
        new{T}(p, s, siglead(p, s))
end

GBElements.cost(g :: SigPoly{T}) where {T} = GBElements.cost(g.polynomial)
GBElements.fullform(g :: SigPoly{T}) where {T} = GBElements.fullform(g.polynomial)

has_signature(:: Type{<: SigPoly}) = true
signature(g :: SigPoly{T}) where {T} = g.signature

function Base.show(
    io :: IO,
    g :: SigPoly{T}
) where {T <: GBElement}
    print(io, g.polynomial, " signature: ", g.signature)
end

#
# Implementation of the AbstractVector interface for SigPolys, based on
# GradedBinomial
#

function Base.size(
    g :: SigPoly{T}
) :: Tuple where {T <: GBElement}
    return size(g.polynomial)
end

function Base.getindex(
    g :: SigPoly{T},
    i :: Int
) :: Int where {T <: GBElement}
    return g.polynomial[i]
end

function Base.setindex(
    g :: SigPoly{T},
    v :: Int,
    i :: Int
) where {T <: GBElement}
    g.polynomial[i] = v
end

function Base.length(
    g :: SigPoly{T}
) :: Int where {T <: GBElement}
    return length(g.polynomial)
end

#
# Dealing with ModuleMonomialOrderings and how to compare signatures
#

"""
An enum containing all types of module monomial orders available for this
signature algorithm implementation.
"""
@enum ModuleMonomialOrder pot_order ltpot_order top_order

"""
Represents a module monomial ordering. This has to store a reference to the
current basis in order to implement an order like lt-pot (Schreyer's order).

This is mutable because, for convenience, it may sometimes be necessary to
change `monomial_order`.
"""
mutable struct ModuleMonomialOrdering{T <: GBElement} <: GBOrder
    monomial_order :: MonomialOrder
    module_order :: ModuleMonomialOrder
    generators :: Vector{SigPoly{T}}

    function ModuleMonomialOrdering(
        costs :: Array{Int, 2},
        module_order :: ModuleMonomialOrder,
        generators :: Vector{SigPoly{T}}
    ) where {T <: GBElement}
        monomial_order = MonomialOrder(costs)
        new{T}(monomial_order, module_order, generators)
    end
end

function ModuleMonomialOrdering(
    costs :: Array{Int, 2},
    module_order :: Symbol,
    generators :: Vector{SigPoly{T}}
) where {T <: GBElement}
    monomial_order = MonomialOrder(costs)
    if module_order == :pot
        mod = pot_order
    elseif module_order == :top
        mod = top_order
    elseif module_order == :ltpot
        mod = ltpot_order
    else
        error("Unknown module monomial order. It must be one of :pot, :top or :ltpot.")
    end
    return ModuleMonomialOrdering(costs, mod, generators)
end

function Orders.is_inverted(
    o :: ModuleMonomialOrdering{T},
    g :: T
) where {T <: GBElement}
    return is_inverted(o.monomial_order, g, cost(g))
end

function Orders.change_ordering!(
    order :: ModuleMonomialOrdering{T},
    new_order :: Array{Int, 2}
) where {T <: GBElement}
    Orders.change_ordering!(order.monomial_order, new_order)
end

"""
Compares two SigPolys as polynomials wrt to the monomial order of `o`.
Ignores signatures and sig-leads.
"""
function Base.lt(
    o :: ModuleMonomialOrdering{T},
    g :: SigPoly{T},
    h :: SigPoly{T}
) :: Bool where {T <: GBElement}
    return Base.lt(o.monomial_order, g.polynomial, h.polynomial)
end

"""
Compares two signatures with respect to the given ModuleMonomialOrdering. Returns
true iff `s1` is strictly smaller than `s2` with respect to `o`.

TODO I wonder if I shouldn't also store the weighted versions of signature.monomial
This could make all orders more efficient, the obvious tradeoff being memory.
"""
function Base.lt(
    o :: ModuleMonomialOrdering{T},
    s1 :: Signature,
    s2 :: Signature
) :: Bool where {T <: GBElement}
    return signature_lt(s1, s2, o)
end

"""
Compares `s1` and `s2` in the position-over-term (pot) order, that is, compare
signatures first by module index and break ties by the given monomial order on the
coefficients.
"""
function pot_compare(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder
) :: Symbol
    if s1.index < s2.index
        return :lt
    elseif s1.index == s2.index
        #c = cmp(monomial_order * s1.monomial, monomial_order * s2.monomial)
        c = cmp(monomial_order, s1.monomial, s2.monomial)
        if c == -1
            return :lt
        elseif c == 0
            return :eq
        else
            return :gt
        end
    end
    return :gt
end

"""
Returns true iff s1 < s2 in the top order.
"""
function pot_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder
) :: Bool
    return pot_compare(s1, s2, monomial_order) == :lt
end

"""
Compares `s1` and `s2` in the term-over-position (top) order, that is, compare
signature first by monomial order on the coefficients breaking ties by the module
index.

Returns :lt, :eq or :gt if `s1` is smaller, equal or greater than `s2`
respectively.
"""
function top_compare(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder
) :: Symbol
    #c = cmp(monomial_order * s1.monomial, monomial_order * s2.monomial)
    c = cmp(monomial_order, s1.monomial, s2.monomial)
    if c == -1
        return :lt
    elseif c == 0
        if s1.index < s2.index
            return :lt
        elseif s1.index == s2.index
            return :eq
        else
            return :gt
        end
    end
    return :gt
end

"""
Returns true iff s1 < s2 in the top order.
"""
function top_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder
) :: Bool
    return top_compare(s1, s2, monomial_order) == :lt
end

"""
Compares `s1` and `s2` in Schreyer's order (aka ltpot). This means first
comparing the leading terms of the images of the signatures in the polynomial ring
and then breaking ties by the pot order.

Returns :lt, :eq or :gt if `s1` is smaller, equal or greater than `s2`
respectively.

TODO try to implement this more efficiently. The way I'm doing this here is very
slow, due to all the calls to image_leading_term...
"""
function ltpot_compare(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder,
    generators :: Vector{SigPoly{T}}
) :: Symbol where {T <: GBElement}
    lt_s1 = image_leading_term(s1, generators)
    lt_s2 = image_leading_term(s2, generators)
    #c = cmp(monomial_order * lt_s1, monomial_order * lt_s2)
    c = cmp(monomial_order, lt_s1, lt_s2)
    if c == -1
        return :lt
    elseif c == 0
        return pot_compare(s1, s2, monomial_order)
    end
    return :gt
end

"""
Returns true iff s1 < s2 in the ltpot (Schreyer) order.
"""
function ltpot_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: MonomialOrder,
    generators :: Vector{SigPoly{T}}
) :: Bool where {T <: GBElement}
    return ltpot_compare(s1, s2, monomial_order, generators) == :lt
end

"""
Returns the leading term image of `s` by the function mapping the module basis
element of index i to the i-th element of generators.
"""
function image_leading_term(
    s :: Signature,
    generators :: Vector{SigPoly{T}}
) :: Vector{Int} where {T <: GBElement}
    @assert s.index <= length(generators)
    gen = generators[s.index].polynomial
    @assert length(s.monomial) == length(gen)
    image_lt = copy(s.monomial)
    gen_lt = leading_term(gen)
    image_lt += gen_lt
    return image_lt
end

"""
Returns :lt, :eq or :gt when `sig1` is respectively smaller, equal or greater than
`sig2` with respect to `module_order`.
"""
function signature_compare(
    sig1 :: Signature,
    sig2 :: Signature,
    module_order :: ModuleMonomialOrdering{T}
) :: Symbol where {T <: GBElement}
    if module_order.module_order == pot_order
        return pot_compare(sig1, sig2, module_order.monomial_order)
    elseif module_order.module_order == top_order
        return top_compare(sig1, sig2, module_order.monomial_order)
    end
    return ltpot_compare(sig1, sig2, module_order.monomial_order,
                         module_order.generators)
end

"""
Returns `sig1` < `sig2` according to the given `module_order`
"""
function signature_lt(
    sig1 :: Signature,
    sig2 :: Signature,
    module_order :: ModuleMonomialOrdering{T}
) :: Bool where {T <: GBElement}
    return signature_compare(sig1, sig2, module_order) == :lt
end

#
# Implementation of CriticalPairs with signatures
#

"""
An S-pair represented sparsely, before building it from binomials explicitly.
Includes a signature to allow the signature-based algorithm to proceed by
increasing signature of S-pairs.
"""
struct SignaturePair <: CriticalPair
    i :: Int
    j :: Int
    signature :: Signature
end

GBElements.first(pair :: SignaturePair) = pair.i
GBElements.second(pair :: SignaturePair) = pair.j

#This is necessary because these SPairs are queued.
function Base.lt(
    o :: ModuleMonomialOrdering{T},
    s1 :: SignaturePair,
    s2 :: SignaturePair
) :: Bool where {T <: GBElement}
    return signature_lt(s1.signature, s2.signature, o)
end

"""
Builds u - v with signature given by `pair`.
"""
function GBElements.build(
    u :: SigPoly{T},
    v :: SigPoly{T},
    pair :: SignaturePair
) :: SigPoly{T} where {T <: GBElement}
    element = u.polynomial - v.polynomial
    signature = pair.signature
    return SigPoly{T}(element, signature)
end

#
# Implementation of a signature basis, useful for passing around in the
# reduction process. It takes the module ordering along in a convenient way.
#

const SigBasis{T} = BinomialSet{SigPoly{T}, ModuleMonomialOrdering{T}}

"""
Computes the signature of the i-branch of the S-pair (i, j). Just call with i and
j inverted to compute the same thing for j.
"""
function spair_signature(
    i :: Int,
    j :: Int,
    gb :: SigBasis{T}
) :: Signature where {T <: GBElement}
    n = length(gb[i])
    coef = Vector{Int}(undef, n)
    for k in 1:n
        lcm = max(gb[i][k], gb[j][k], 0)
        #For efficiency, we can already add the signature of i here.
        coef[k] = lcm - max(0, gb[i][k]) + gb[i].signature.monomial[k]
    end
    return Signature(gb[i].signature.index, coef)
end

"""
Creates an SPair S(i, j) if it is regular, otherwise returns `nothing`.
"""
function regular_spair(
    i :: Int,
    j :: Int,
    gb :: SigBasis{T};
    comparator :: Union{Comparator{SigLead}, Nothing} = nothing
) :: Union{SignaturePair, Nothing} where {T <: GBElement}
    #Compare sig-leads to determine the signature efficiently
    if !isnothing(comparator)
        comp = compare(comparator, i, j)
    else
        comp = signature_compare(gb[i].siglead, gb[j].siglead, order(gb))
    end
    if comp == :lt #Signature is based on the j-branch
        signature = spair_signature(j, i, gb)
        return SignaturePair(i, j, signature)
    elseif comp == :gt #Signature is based on the i-branch
        signature = spair_signature(i, j, gb)
        return SignaturePair(i, j, signature)
    end
    #This S-pair is singular, eliminate it
    return nothing
end

#
# Implementation of the GBElement interface for SigPolys
#

GBElements.is_zero(g :: SigPoly{GBElement}) = GBElements.is_zero(g.polynomial)
GBElements.head(g :: SigPoly{GBElement}) = head(g.polynomial)
GBElements.supports(g :: SigPoly{GBElement}) = GBElements.supports(g.polynomial)

"""
Checks whether g is singular top-reducible by some reducer with signature
reducer_sig. Assumes this reducer divides g.

TODO singular_top_reducible and signature_reducible both compare the same
signatures. I could call the same function signature_compare in both!
In fact, given the situation I use this, I should just call signature_compare
there and be done with it!
"""
function GBElements.singular_top_reducible(
    g :: SigPoly{T},
    reducer_sig :: Signature
) :: Bool where {T <: GBElement}
    return isequal(g.signature, reducer_sig)
end

"""
Checks whether g is signature reducible by a reducer with reducer_sig. This means
the reducer (after multiplying by a lt quotient) has lower signature than g wrt
the module monomial ordering.
"""
function GBElements.signature_reducible(
    g :: SigPoly{T},
    reducer_sig :: Signature,
    gb :: SigBasis{T}
) :: Bool where {T <: GBElement}
    return Base.lt(order(gb), reducer_sig, g.signature)
end

"""
Computes g -= h in-place, assumes this is a regular reduction, so that
sig(g - h) = sig(g).

Returns true iff `g` reduced to zero.
"""
function GBElements.reduce!(
    g :: SigPoly{T},
    h :: SigPoly{T},
    order :: GBOrder;
    negative :: Bool = false
) :: Bool where {T <: GBElement}
    #No need to update signature, as we assume it won't change
    reduced_to_zero = GBElements.reduce!(g.polynomial, h.polynomial, order, negative=negative)
    #TODO this is necessary, but is inefficient as is. How can I do this cleanly and efficiently?
    #TODO we ignore the negative parameter here. Currently, it is only necessary for reduced
    #basis, and this is only done at the end of the computation. At that point, sigleads
    #are irrelevant.
    for i in 1:length(g.siglead.monomial)
        lt_i = max(g[i], 0)
        g.siglead.monomial[i] = g.signature.monomial[i] - lt_i;
    end
    return reduced_to_zero
end

"""
Checks a reduction condition for implicit variables in the case
T == GradedBinomial. Currently unused, because GradedBinomials are not supported
with signatures.
"""
function GBElements.degree_reducible(
    g :: SigPoly{T},
    h :: SigPoly{T};
    negative :: Bool = false
) :: Bool where {T <: GBElement}
    return GBElements.degree_reducible(
        g.polynomial, h.polynomial, negative=negative
    )
end

"""
Returns the Koszul signature of the pair of gb indexed by (i, j), if it is
regular. Otherwise, returns nothing.
"""
function koszul(
    i :: Int,
    j :: Int,
    gb :: SigBasis{T};
    comparator :: Union{Comparator{SigLead}, Nothing} = nothing
) :: Union{Signature, Nothing} where {T <: GBElement}
    #Compare sig-leads to find out the largest signature
    if !isnothing(comparator)
        comp = compare(comparator, i, j)
    else
        comp = signature_compare(gb[i].siglead, gb[j].siglead, order(gb))
    end
    if comp == :lt
        return leading_term(gb[i]) * gb[j].signature
    elseif comp == :gt
        return leading_term(gb[j]) * gb[i].signature
    end
    #This Koszul syzygy is singular, ignore it
    return nothing
end

end
