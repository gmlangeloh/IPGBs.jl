"""
Implementation of the basic signature-polynomial pairs structures.

TODO implement a getfilter method for SigPolys
"""
module SignaturePolynomials
export Signature, SigPoly, SigBasis, ModuleMonomialOrdering, SPair,
    ModuleMonomialOrder, regular_spair, build_spair, projection, is_zero,
    divides

using IPGBs.GBElements
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.SupportTrees

"""
An enum containing all types of module monomial orders available for my
signature algorithm implementation.
"""
@enum ModuleMonomialOrder begin
    pot = -1
    ltpot = 0
    top = 1
end

"""
A signature, or module monomial, is represented by a monomial with an index,
where the monomial is the coefficient of the module basis vector of the given
index.
"""
struct Signature
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
    if all(i -> s1.monomial[i] <= s2.monomial[i], 1:n)
        return true
    end
    return false
end

function Base.:*(
    monomial :: Vector{Int},
    s :: Signature
) :: Signature
    new_monomial = zeros(Int, length(monomial))
    @assert length(new_monomial) == length(s.monomial)
    for i in 1:length(monomial)
        new_monomial[i] = monomial[i] + s.monomial[i]
    end
    return Signature(s.index, new_monomial)
end

"""
Represents a polynomial with its signature. Used in signature-based algorithms.
"""
struct SigPoly{T <: GBElement} <: GBElement
    polynomial :: T
    signature :: Signature
end

function Base.show(
    io :: IO,
    g :: SigPoly{T}
) where {T <: GBElement}
    print(io, g.polynomial, " signature: ", g.signature)
end

"""
Returns the representation of this SigPoly as a vector of integers, dropping its
signature.
"""
function projection(
    sigpoly :: SigPoly{T}
) :: Vector{Int} where {T <: GBElement}
    return sigpoly.polynomial.element
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

struct ModuleMonomialOrdering{T <: GBElement} <: Base.Ordering
    monomial_order :: Array{Int, 2}
    module_order :: ModuleMonomialOrder
    #TODO I should make sure this is a reference to the same GB object that
    #is incrementally built in my signature algorithm
    generators :: Vector{SigPoly{T}}
end

function Base.lt(
    o :: ModuleMonomialOrdering,
    s1 :: Signature,
    s2 :: Signature
) :: Bool
    return signature_lt(s1, s2, o.monomial_order, o.generators, o.module_order)
end

"""
Returns true iff s1 < s2 in the position-over-term (pot) order, that is,
compare signatures first by module index and break ties by the given
monomial order on the coefficients.
"""
function pot_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: Array{Int, 2}
) :: Bool
    if s1.index < s2.index
        return true
    elseif s1.index == s2.index
        return monomial_order * s1.monomial < monomial_order * s2.monomial
    end
    return false
end

"""
Returns true iff s1 < s2 in the term-over-position (top) order, that is,
compare signature first by monomial order on the coefficients breaking
ties by the module index.
"""
function top_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: Array{Int, 2}
) :: Bool
    weighted_s1 = monomial_order * s1.monomial
    weighted_s2 = monomial_order * s2.monomial
    if weighted_s1 < weighted_s2
        return true
    elseif weighted_s1 == weighted_s2
        return s1.index < s2.index
    end
    return false
end

"""
Returns true iff s1 < s2 in Schreyer's order (aka ltpot). This means
first comparing the leading terms of the images of the signatures in
the polynomial ring and then breaking ties by the pot order.
"""
function ltpot_lt(
    s1 :: Signature,
    s2 :: Signature,
    monomial_order :: Array{Int, 2},
    generators :: Vector{SigPoly{T}}
) :: Bool where {T <: GBElement}
    lt_s1 = image_leading_term(s1, generators)
    lt_s2 = image_leading_term(s2, generators)
    weighted_lts1 = monomial_order * lt_s1
    weighted_lts2 = monomial_order * lt_s2
    if weighted_lts1 < weighted_lts2
        return true
    elseif weighted_lts1 == weighted_lts2
        return pot_lt(s1, s2, monomial_order)
    end
    return false
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
Returns sig1 < sig2 according to the given module_order and monomial_order.
"""
function signature_lt(
    sig1 :: Signature,
    sig2 :: Signature,
    monomial_order :: Array{Int, 2},
    generators :: Vector{SigPoly{T}},
    module_order :: ModuleMonomialOrder
) :: Bool where {T <: GBElement}
    if module_order == pot
        return pot_lt(sig1, sig2, monomial_order)
    elseif module_order == ltpot
        return ltpot_lt(sig1, sig2, monomial_order, generators)
    end
    #module_order == top
    return top_lt(sig1, sig2, monomial_order)
end

"""
An S-pair represented sparsely, before building it from binomials explicitly.
Includes a signature to allow the signature-based algorithm to proceed by
increasing signature of S-pairs.
"""
struct SPair
    i :: Int
    j :: Int
    signature :: Signature
end

function Base.lt(
    o :: ModuleMonomialOrdering,
    s1 :: SPair,
    s2 :: SPair
) :: Bool
    return signature_lt(s1.signature, s2.signature, o.monomial_order,
                        o.generators, o.module_order)
end

#
# Implementation of a signature basis, useful for passing around in the
# reduction process. It takes the module ordering along in a convenient way.
#

struct SigBasis{T <: GBElement} <: AbstractVector{SigPoly{T}}
    basis :: Vector{SigPoly{T}}
    module_ordering :: ModuleMonomialOrdering
    reduction_tree :: SupportTree{SigPoly{T}}
end

function Base.size(
    gb :: SigBasis{T}
) :: Tuple where {T <: GBElement}
    return size(gb.basis)
end

function Base.getindex(
    gb :: SigBasis{T},
    i :: Int
) :: SigPoly{T} where {T <: GBElement}
    return gb.basis[i]
end

function Base.setindex(
    gb :: SigBasis{T},
    g :: SigPoly{T},
    i :: Int
) where {T <: GBElement}
    gb.basis[i] = g
end

function Base.length(
    gb :: SigBasis{T}
) :: Int where {T <: GBElement}
    return length(gb.basis)
end

function Base.push!(
    gb :: SigBasis{T},
    g :: SigPoly{T}
) where {T <: GBElement}
    push!(gb.basis, g)
    push!(gb.module_ordering.generators, g)
    addbinomial!(gb.reduction_tree, gb[length(gb)])
end

#
# Implementation of the GBElement interface for SigPolys
#

GBElements.is_zero(g :: SigPoly{GBElement}) = GBElements.is_zero(g.polynomial)

"""
Checks whether a * reducer.signature < g.signature wrt a given module monomial
order, where a is a monomial such that a * lt(reducer) == lt(g).

TODO I still need to pass on the fact that the thing is singular s-reducible
which will happen in the SupportTree somehow...
"""
function GBElements.regular_reducible(
    reducer :: SigPoly{T},
    g :: SigPoly{T},
    gb :: SigBasis{T}
) :: Bool where {T <: GBElement}
    #Compute the multiplication factor
    n = length(reducer)
    factor = zeros(Int, n)
    for i in head(g.polynomial)
        factor[i] = g[i] - reducer[i]
        @assert factor[i] >= 0
    end

    @show factor * reducer.signature
    @show g.signature

    #DEBUGGING
    bugged = [
        [-1, 1, 0, 0, 0, 0],
        [1, 0, 0, -1, 1, 0],
        [0, -1, 1, 0, 1, -1]
    ]
    if g.polynomial.element in bugged
        println("bug case?")
        @show reducer reducer.signature
        @show g.signature
        @show factor * reducer.signature
    end

    return Base.lt(gb.module_ordering, factor * reducer.signature,
                   g.signature)
end

"""
Computes g -= h in-place, assumes this is a regular reduction, so that
sig(g - h) = sig(g).

Returns true iff `g` reduced to zero.
"""
function GBElements.reduce!(
    g :: SigPoly{T},
    h :: SigPoly{T}
) :: Bool where {T <: GBElement}
    #No need to update signature, as we assume it won't change
    return GBElements.reduce!(g.polynomial, h.polynomial)
end

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
Builds concrete S-pair from an `SPair` struct.
"""
function build_spair(
    spair :: SPair,
    generators :: SigBasis{T}
) :: SigPoly{T} where {T <: GBElement}
    g_i = generators[spair.i].polynomial
    g_j = generators[spair.j].polynomial
    if g_i.cost < g_j.cost
        s = g_j - g_i
    elseif g_i.cost > g_j.cost
        s = g_i - g_j
    else
        #Do a tiebreaker
        #TODO Maybe I should pass the whole matrix C here for tiebreaking...
        if GBElements.lt_tiebreaker(g_i, g_j)
            s = g_j - g_i
        else
            s = g_i - g_j
        end
    end
    return SigPoly(s, spair.signature)
end

"""
Returns (a, b) for SPair(i, j) = ag_i + bg_j.
"""
function spair_coefs(
    i :: Int,
    j :: Int,
    gb :: SigBasis{T}
) :: Tuple{Vector{Int}, Vector{Int}} where {T <: GBElement}
    n = length(gb[i].polynomial)
    i_coef = zeros(Int, n)
    j_coef = zeros(Int, n)
    for k in 1:n
        lcm = max(gb[i].polynomial[k], gb[j].polynomial[k], 0)
        i_coef[k] = lcm - max(0, gb[i].polynomial[k])
        j_coef[k] = lcm - max(0, gb[j].polynomial[k])
    end
    return i_coef, j_coef
end

"""
Creates an SPair S(i, j) if it is regular, otherwise returns `nothing`.
"""
function regular_spair(
    i :: Int,
    j :: Int,
    gb :: SigBasis{T}
) :: Union{SPair, Nothing} where {T <: GBElement}
    i_coef, j_coef = spair_coefs(i, j, gb)
    i_sig = i_coef * gb[i].signature
    j_sig = j_coef * gb[j].signature
    if i_sig == j_sig #S-pair is singular, eliminate
        return nothing
    end #otherwise s-pair is regular, generate it
    sig_lt = signature_lt(i_sig, j_sig, gb.module_ordering.monomial_order,
                          gb.basis, gb.module_ordering.module_order)
    if sig_lt
        sig = j_sig
    else
        sig = i_sig
    end
    return SPair(i, j, sig)
end

end
