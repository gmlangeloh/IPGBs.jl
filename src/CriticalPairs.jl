module CriticalPairs

export CriticalPair, BinomialPair, SignaturePair, first, second

using IPGBs.GBElements
using IPGBs.SignaturePolynomials

abstract type CriticalPair end

#CriticalPair interface
first(:: CriticalPair) :: Int = error("Not implemented.")
second(:: CriticalPair) :: Int = error("Not implemented.")

#
# Basic CriticalPair, used in Buchberger's algorithm
#

struct BinomialPair <: CriticalPair
    i :: Int
    j :: Int
end

first(pair :: BinomialPair) = pair.i
second(pair :: BinomialPair) = pair.j

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

first(pair :: SignaturePair) = pair.i
second(pair :: SignaturePair) = pair.j

#This is necessary because these SPairs are queued.
function Base.lt(
    o :: ModuleMonomialOrdering{T},
    s1 :: SignaturePair,
    s2 :: SignaturePair
) :: Bool where {T <: GBElement}
    return signature_lt(s1.signature, s2.signature, o.monomial_order,
                        o.generators, o.module_order)
end

end
