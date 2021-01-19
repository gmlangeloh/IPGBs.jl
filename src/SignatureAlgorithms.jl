"""
An implementation of a Signature-based algorithm for Gröbner bases of toric
ideals.
"""
module SignatureAlgorithms
export SignatureAlgorithm

using DataStructures

using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.SignaturePolynomials
using IPGBs.SupportTrees

using IPGBs.GBAlgorithms

#
# Syzygy division
#

"""
Stores one SupportTree for each signature index. This enables fast divisor
queries for Signatures.
"""
mutable struct SyzygySet
    syzygies :: Vector{SupportTree{Signature}}
    total_syzygies :: Int
    SyzygySet() = new(SupportTree{Signature}[], 0)
end

function initialize!(
    syzygies :: SyzygySet,
    n :: Int
)
    for i in 1:n
        t = support_tree(Signature[])
        push!(syzygies.syzygies, t)
    end
    syzygies.total_syzygies = n
end

function add_syzygy!(
    syzygies :: SyzygySet,
    new_syzygy :: Signature
)
    addbinomial!(syzygies.syzygies[new_syzygy.index], new_syzygy)
    syzygies.total_syzygies += 1
end

"""
Returns true iff there is a syzygy in `syzygies` that divides `s`. Uses
SupportTrees for efficient divisor search.

TODO use this for the signature criterion
"""
function divided_by(
    s :: Signature,
    syzygies :: SyzygySet
) :: Bool
    search_tree = syzygies.syzygies[s.index]
    #I need to pass a signature array here because the function below uses the
    #metadata on a BinomialSet, when it is available. An empty vector with no
    #additional data works fine, though.
    divisor = find_reducer(s, Signature[], search_tree)
    return !isnothing(divisor)
end

#
# Definition of Signature Algorithms themselves.
#

#Type alias for my SignaturePair heaps
#TODO maybe I should use some other name though, due to the Koszul heaps
const SigHeap{T} = BinaryHeap{SignaturePair, ModuleMonomialOrdering{T}};

mutable struct SignatureAlgorithm{T} <: GBAlgorithm
    basis :: SigBasis{T}
    heap :: SigHeap{T}
    #syzygies :: Vector{Signature}
    syzygies :: SyzygySet
    previous_signature :: Union{Signature, Nothing}
    stats :: GBStats

    function SignatureAlgorithm(T :: Type, C :: Array{Int, 2})
        syzygies = SyzygySet()
        generators = SigPoly{T}[]
        #For now, fix the Schreyer order as module monomial order
        order = ModuleMonomialOrdering(C, SignaturePolynomials.ltpot, generators)
        basis = BinomialSet(generators, order)
        heap = BinaryHeap{SignaturePair}(order, [])
        #TODO create additional stats fields for Signatures
        #This should include counting the number of pairs eliminated by each
        #criterion
        stats = GBStats()
        stats.stats["eliminated_by_signature"] = 0
        stats.stats["eliminated_by_gcd"] = 0
        stats.stats["eliminated_by_duplicate"] = 0
        new{T}(basis, heap, syzygies, nothing, stats)
    end
end

GBAlgorithms.stats(algorithm :: SignatureAlgorithm) = algorithm.stats
GBAlgorithms.data(algorithm :: SignatureAlgorithm) = algorithm.stats.stats
GBAlgorithms.current_basis(algorithm :: SignatureAlgorithm) = algorithm.basis

function GBAlgorithms.next_pair!(
    algorithm :: SignatureAlgorithm{T}
) :: Union{SignaturePair, Nothing} where {T <: GBElement}
    if isempty(algorithm.heap)
        return nothing
    end
    return pop!(algorithm.heap)
end

function early_pair_elimination(
    algorithm :: SignatureAlgorithm{T},
    pair :: SignaturePair
) :: Bool where {T <: GBElement}
    #GCD criterion
    if is_support_reducible(GBElements.first(pair), GBElements.second(pair),
                            current_basis(algorithm))
        data(algorithm)["eliminated_by_gcd"] += 1
        return true
    end
    return false
end

function GBAlgorithms.late_pair_elimination(
    algorithm :: SignatureAlgorithm{T},
    pair :: SignaturePair
) :: Bool where {T <: GBElement}
    #Skip repeated signatures
    if !isnothing(algorithm.previous_signature) &&
        isequal(algorithm.previous_signature, pair.signature)
        data(algorithm)["eliminated_by_duplicate"] += 1
        return true
    end
    algorithm.previous_signature = pair.signature
    #Signature criterion
    if signature_criterion(pair, algorithm.syzygies)
        data(algorithm)["eliminated_by_signature"] += 1
        return true
    end
    #GCD criterion
    if is_support_reducible(GBElements.first(pair), GBElements.second(pair),
                            current_basis(algorithm))
        data(algorithm)["eliminated_by_gcd"] += 1
        return true
    end
    return false
end

function GBAlgorithms.process_zero_reduction!(
    algorithm :: SignatureAlgorithm{T},
    syzygy_element :: SigPoly{T}
) where {T <: GBElement}
    add_syzygy!(algorithm.syzygies, signature(syzygy_element))
    #push!(algorithm.syzygies, signature(syzygy_element))
    data(algorithm)["zero_reductions"] += 1
end

function GBAlgorithms.update!(
    algorithm :: SignatureAlgorithm{T},
    g :: SigPoly{T}
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    update_queue!(algorithm)
    update_syzygies!(algorithm)
    data(algorithm)["max_basis_size"] = max(data(algorithm)["max_basis_size"],
                                            length(current_basis(algorithm)))
end

"""
Adds new SignaturePairs to algorithm.heap. Applies early elimination of S-pairs
when possible.
"""
function update_queue!(
    algorithm :: SignatureAlgorithm{T}
) where {T <: GBElement}
    gb = current_basis(algorithm)
    n = length(gb)
    for i in 1:(n-1)
        sp = regular_spair(i, n, gb)
        if !isnothing(sp) && !early_pair_elimination(algorithm, sp)
            push!(algorithm.heap, sp)
            data(algorithm)["queued_pairs"] += 1
        end
    end
end

"""
Adds the Koszul syzygies corresponding to the newest element of gb to the
syzygy list.
"""
function update_syzygies!(
    algorithm :: SignatureAlgorithm{T}
) where {T <: GBElement}
    gb = current_basis(algorithm)
    n = length(gb)
    for i in 1:(n-1)
        add_syzygy!(algorithm.syzygies, koszul(i, n, gb))
        #push!(algorithm.syzygies, koszul(i, n, gb))
    end
end

function change_ordering!(
    algorithm :: SignatureAlgorithm{T},
    new_monomial_order :: Array{Int, 2}
) where {T <: GBElement}
    SignaturePolynomials.change_ordering!(
        algorithm.basis.order, new_monomial_order
    )
end

#TODO almost all of this logic is still duplicated in Buchberger.jl. Can I do
#this better?
function GBAlgorithms.initialize!(
    algorithm :: SignatureAlgorithm{T},
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    # This assumes binary constraints
    change_ordering!(algorithm, C)
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    #Need to initialize the SyzygySet before adding new elements in case we
    #already start adding Koszul syzygies
    initialize!(algorithm.syzygies, num_gens)
    num_vars = size(A, 2)
    coef = zeros(Int, num_vars) #Coefficient of the signatures of these generators
    j = 1
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u)
        if !isnothing(e)
            s = Signature(j, coef)
            GBAlgorithms.update!(algorithm, SigPoly{T}(e, s))
            j += 1
        end
    end
end

#
# Signature S-pair elimination criteria
#

"""
Returns true iff spair is eliminated by the signature criterion, i.e. there
is a known syzygy that divides its signature.
"""
function signature_criterion(
    spair :: SignaturePair,
    #syzygies :: Vector{Signature}
    syzygies :: SyzygySet
) :: Bool
    #for syz in syzygies
    #    if divides(syz, spair.signature)
    #        return true
    #    end
    #end
    #return false
    return divided_by(spair.signature, syzygies)
end

end
