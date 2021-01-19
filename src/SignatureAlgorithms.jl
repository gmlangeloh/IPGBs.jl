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
const SigHeap{T} = BinaryHeap{SignaturePair, ModuleMonomialOrdering{T}}
const KoszulHeap{T} = BinaryHeap{Signature, ModuleMonomialOrdering{T}}

mutable struct SignatureAlgorithm{T} <: GBAlgorithm
    basis :: SigBasis{T}
    heap :: SigHeap{T}
    koszul :: KoszulHeap{T}
    syzygies :: SyzygySet
    previous_signature :: Union{Signature, Nothing}
    stats :: GBStats

    function SignatureAlgorithm(T :: Type, C :: Array{Int, 2}, mod_order :: Symbol)
        syzygies = SyzygySet()
        generators = SigPoly{T}[]
        #For now, fix the Schreyer order as module monomial order
        if mod_order == :pot
            module_order = pot_order
        elseif mod_order == :ltpot
            module_order = ltpot_order
        else
            module_order = top_order
        end
        order = ModuleMonomialOrdering(C, module_order, generators)
        basis = BinomialSet(generators, order)
        heap = BinaryHeap{SignaturePair}(order, [])
        koszul = BinaryHeap{Signature}(order, [])
        stats = GBStats()
        stats.stats["eliminated_by_early_signature"] = 0
        stats.stats["eliminated_by_late_signature"] = 0
        stats.stats["eliminated_by_gcd"] = 0
        stats.stats["eliminated_by_duplicate"] = 0
        stats.stats["eliminated_by_early_koszul"] = 0
        stats.stats["eliminated_by_late_koszul"] = 0
        new{T}(basis, heap, koszul, syzygies, nothing, stats)
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
    #Signature criterion
    if signature_criterion(pair, algorithm.syzygies)
        data(algorithm)["eliminated_by_early_signature"] += 1
        return true
    end
    if koszul_criterion(pair, algorithm.koszul, algorithm.basis.order)
        data(algorithm)["eliminated_by_early_koszul"] += 1
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
        data(algorithm)["eliminated_by_late_signature"] += 1
        return true
    end
    if koszul_criterion(pair, algorithm.koszul, algorithm.basis.order)
        data(algorithm)["eliminated_by_late_koszul"] += 1
        return true
    end
    #GCD criterion
    #Unnecessary to apply GCD criterion both early and late. Just use one of them.
    #if is_support_reducible(GBElements.first(pair), GBElements.second(pair),
    #                        current_basis(algorithm))
    #    data(algorithm)["eliminated_by_late_gcd"] += 1
    #    return true
    #end
    return false
end

function GBAlgorithms.process_zero_reduction!(
    algorithm :: SignatureAlgorithm{T},
    syzygy_element :: SigPoly{T}
) where {T <: GBElement}
    add_syzygy!(algorithm.syzygies, signature(syzygy_element))
    data(algorithm)["zero_reductions"] += 1
end

function GBAlgorithms.update!(
    algorithm :: SignatureAlgorithm{T},
    g :: SigPoly{T},
    pair :: Union{SignaturePair, Nothing} = nothing
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    update_queue!(algorithm)
    #update_syzygies!(algorithm)
    if !isnothing(pair)
        push!(algorithm.koszul, koszul(
            GBElements.first(pair), GBElements.second(pair),
            current_basis(algorithm)))
    end
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
        #add_syzygy!(algorithm.syzygies, koszul(i, n, gb))
        push!(algorithm.koszul, koszul(i, n, gb))
    end
end

"""
Changes to a new monomial ordering. This may be necessary in the beginning of
the run in case the monomial order has suffered some transformation (e.g. new
variables were added due to normalization) and thus these changes must be also
done in the algorithm's internal order.
"""
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
    syzygies :: SyzygySet
) :: Bool
    return divided_by(spair.signature, syzygies)
end

#TODO this isn't eliminating pretty much anything at all! Investigate
function koszul_criterion(
    spair :: SignaturePair,
    koszul :: KoszulHeap{T},
    order :: ModuleMonomialOrdering{T}
) :: Bool where {T <: GBElement}
    if isempty(koszul) #No known Koszul signatures, skip criterion
        return false
    end
    #Pop all Koszul signatures that are smaller than the current spair signature
    #Because signatures are processed in increasing order, if these signatures
    #are useless now, they will always be useless later
    k_sig = Base.first(koszul)
    while Base.lt(order, k_sig, spair.signature)
        pop!(koszul)
        if isempty(koszul)
            return false
        end
        k_sig = Base.first(koszul)
    end
    #If current Koszul signature and spair signature coincide, the spair may be
    #eliminated
    if Base.isequal(k_sig, spair.signature)
        return true
    end
    return false
end

end
