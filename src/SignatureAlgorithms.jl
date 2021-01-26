"""
An implementation of a Signature-based algorithm for Gröbner bases of toric ideals.
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
using IPGBs.TriangleHeaps

using IPGBs.GBAlgorithms

#
# Syzygy division
#

"""
Stores one SupportTree for each signature index. This enables fast divisor
queries for Signatures.
"""
mutable struct SignatureSet
    signatures :: Vector{SupportTree{Signature}}
    total_signatures :: Int
    SignatureSet() = new(SupportTree{Signature}[], 0)
end

function initialize!(
    signatures :: SignatureSet,
    n :: Int
)
    for i in 1:n
        t = support_tree(Signature[])
        push!(signatures.signatures, t)
    end
    signatures.total_signatures = n
end

function add_signature!(
    signatures :: SignatureSet,
    new_signature :: Signature
)
    addbinomial!(signatures.signatures[new_signature.index], new_signature)
    signatures.total_signatures += 1
end

"""
Returns true iff there is a signature in `signatures` that divides `s`. Uses
SupportTrees for efficient divisor search.
"""
function divided_by(
    s :: Signature,
    signatures :: SignatureSet
) :: Bool
    search_tree = signatures.signatures[s.index]
    #I need to pass a signature array here because the function below uses the
    #metadata on a BinomialSet, when it is available. An empty vector with no
    #additional data works fine, though.
    divisor = find_reducer(s, Signature[], search_tree)
    return !isnothing(divisor)
end

#
# Definition of Signature Algorithms themselves.
#

const KoszulHeap{T} = BinaryHeap{Signature, ModuleMonomialOrdering{T}}

mutable struct SignatureAlgorithm{T} <: GBAlgorithm
    basis :: SigBasis{T}
    heap :: TriangleHeap{T, UInt32}
    koszul :: KoszulHeap{T}
    syzygies :: SignatureSet
    previous_signature :: Union{Signature, Nothing}
    stats :: GBStats

    function SignatureAlgorithm(T :: Type, C :: Array{Int, 2}, mod_order :: Symbol)
        syzygies = SignatureSet()
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
        heap = TriangleHeap{T, UInt32}(basis, order)
        koszul = BinaryHeap{Signature}(order, [])
        #Initialize all statistics collected by a SignatureAlgorithm
        stats = GBStats()
        stats.stats["eliminated_by_early_signature"] = 0
        stats.stats["eliminated_by_late_signature"] = 0
        stats.stats["eliminated_by_gcd"] = 0
        stats.stats["eliminated_by_duplicate"] = 0
        stats.stats["eliminated_by_early_koszul"] = 0
        stats.stats["eliminated_by_late_koszul"] = 0
        stats.stats["eliminated_by_base_divisors"] = 0
        stats.stats["eliminated_by_low_base_divisors"] = 0
        stats.stats["eliminated_by_high_base_divisors"] = 0
        new{T}(basis, heap, koszul, syzygies, nothing, stats)
    end
end

GBAlgorithms.stats(algorithm :: SignatureAlgorithm) = algorithm.stats
GBAlgorithms.data(algorithm :: SignatureAlgorithm) = algorithm.stats.stats
GBAlgorithms.current_basis(algorithm :: SignatureAlgorithm) = algorithm.basis

"""
Generates the next pair to be processed by `algorithm` by picking the first
element of its priority queue. Returns `nothing` if there are no more S-pairs
to be processed and the algorithm may thus terminate.

Pairs are processed in order of increasing signature.
"""
function GBAlgorithms.next_pair!(
    algorithm :: SignatureAlgorithm{T}
) :: Union{SignaturePair, Nothing} where {T <: GBElement}
    if isempty(algorithm.heap)
        return nothing
    end
    return pop!(algorithm.heap)
end

"""
Attempts to eliminate the pair (i, j) before even creating it as a SignaturePair
(thus, without computing its signature).

Returns true iff (i, j) may be eliminated at this point.
"""
function very_early_pair_elimination(
    algorithm :: SignatureAlgorithm{T},
    i :: Int,
    j :: Int
) :: Bool where {T <: GBElement}
    #GCD criterion
    if is_support_reducible(i, j, current_basis(algorithm))
        data(algorithm)["eliminated_by_gcd"] += 1
        return true
    end
    if base_divisors_criterion(i, j, algorithm)
        data(algorithm)["eliminated_by_base_divisors"] += 1
    end
    return false
end

"""
Attempts to eliminate `pair` before putting it into the priority queue, but after
computing its signature.

Returns true iff `pair` may be eliminated at this point.
"""
function early_pair_elimination(
    algorithm :: SignatureAlgorithm{T},
    pair :: SignaturePair
) :: Bool where {T <: GBElement}
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

"""
Attempts to eliminate `pair` before building it as an explicit GBElement, after
retrieving it from the priority queue for processing.

Returns true iff `pair` may be eliminated at this point.
"""
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
    return false
end

"""
Adds `syzygy_element` to the set of known syzygies of `algorithm`.
"""
function GBAlgorithms.process_zero_reduction!(
    algorithm :: SignatureAlgorithm{T},
    syzygy_element :: SigPoly{T}
) where {T <: GBElement}
    add_signature!(algorithm.syzygies, signature(syzygy_element))
    data(algorithm)["zero_reductions"] += 1
end

"""
Adds the element `g` to the current basis maintained by `algorithm`. Uses `pair`,
if available, to generate new Koszul syzygies. `pair` is assumed to be the
SignaturePair that was used to generate `g`.
"""
function GBAlgorithms.update!(
    algorithm :: SignatureAlgorithm{T},
    g :: SigPoly{T},
    pair :: Union{SignaturePair, Nothing} = nothing
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    update_queue!(algorithm)
    #update_syzygies!(algorithm)
    if !isnothing(pair)
        k = koszul(GBElements.first(pair), GBElements.second(pair),
                   current_basis(algorithm))
        if !isnothing(k)
            push!(algorithm.koszul, k)
        end
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
    batch = SignaturePair[]
    for i in 1:(n-1)
        if very_early_pair_elimination(algorithm, i, n)
            continue
        end
        sp = regular_spair(i, n, gb)
        if !isnothing(sp) && !early_pair_elimination(algorithm, sp)
            push!(batch, sp)
            data(algorithm)["queued_pairs"] += 1
        end
    end
    push_batch!(algorithm.heap, batch)
end

"""
Adds the Koszul syzygies corresponding to the newest element of gb to the
syzygy list.

Currently unused, as this is inefficient. It's not necessary to add all Koszul
syzygies.
"""
function update_syzygies!(
    algorithm :: SignatureAlgorithm{T}
) where {T <: GBElement}
    gb = current_basis(algorithm)
    n = length(gb)
    for i in 1:(n-1)
        k = koszul(i, n, gb)
        if !isnothing(k)
            push!(algorithm.koszul, k)
        end
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

"""
Initializes all data in `algorithm`. Must be called before using any of the data
inside for GB computations.

TODO almost all of this logic is still duplicated in Buchberger.jl. Can I do
#this better?
"""
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
    syzygies :: SignatureSet
) :: Bool
    return divided_by(spair.signature, syzygies)
end

"""
Applies the Koszul criterion to `spair`, returning true iff it can be eliminated
by this criterion. The main idea is checking whether there exists a regular Koszul
syzygy with the same signature as `spair`. If this happens, `spair` may be
eliminated.

TODO this isn't eliminating pretty much anything at all! Investigate
"""
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

"""
Applies the Base Divisors Criterion presented in Roune and Stillman (2012).
"""
function base_divisors_criterion(
    i :: Int,
    j :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}
    if high_base_divisors_criterion(i, j, algorithm)
        data(algorithm)["eliminated_by_high_base_divisors"] += 1
        return true
    end
    const num_low_divisors = 2
    if low_base_divisors_criterion(i, j, algorithm)
        data(algorithm)["eliminated_by_low_base_divisors"] += 1
        return true
    end
    return false
end

"""
The Base divisor criterion for high-ratio elements. An element a is high-ratio
if lt(a) | lt(b), where b is the new element being added.
"""
function high_base_divisors_criterion(
    i :: Int,
    j :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}

end

"""
The Base divisor criterion for low-ratio elements. An element a is low-ratio if
signature(a) | signature(b), where b is the new element being added.
"""
function low_base_divisors_criterion(
    i :: Int,
    j :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}

end

end
