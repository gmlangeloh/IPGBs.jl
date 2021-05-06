"""
An implementation of a Signature-based algorithm for Gröbner bases of toric ideals.
"""
module SignatureAlgorithms
export SignatureAlgorithm

using DataStructures

using IPGBs.FastBitSets
using IPGBs.FastComparator
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.SignaturePolynomials
using IPGBs.Statistics
using IPGBs.SupportTrees
using IPGBs.TriangleHeaps

using IPGBs.GBAlgorithms

#
# Just a simple struct to make the implementation of the Base Divisor criterion simpler
#

struct LowBaseDiv
    index :: Int
    monomial :: Vector{Int}

    function LowBaseDiv(index, a :: SigPoly{T}, b :: SigPoly{T}) where {T}
        n = length(a)
        p = Vector{Int}(undef, n)
        v = Vector{Int}(undef, n)
        lt_a = leading_term(a)
        lt_b = leading_term(b)
        for i in 1:n
            p[i] = lt_a[i] - a.signature.monomial[i] + b.signature.monomial[i]
            if lt_b[i] <= p[i]
                v[i] = typemax(Int)
            else
                v[i] = max(p[i], lt_a[i])
            end
        end
        new(index, v)
    end
end

#
# Signature division (useful for the signature criterion, for example)
#

"""
Stores one SupportTree for each signature index. This enables fast divisor
queries for Signatures.
"""
mutable struct SignatureSet
    signatures :: Vector{SupportTree{Signature}}
    basis_indices :: Dict{Signature, Int}
    total_signatures :: Int
    SignatureSet() = new(SupportTree{Signature}[], Dict{Signature, Int}(), 0)
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
    new_signature :: Signature,
    index :: Union{Int, Nothing} = nothing
)
    addbinomial!(signatures.signatures[new_signature.index], new_signature)
    if !isnothing(index)
        signatures.basis_indices[new_signature] = index
    end
    signatures.total_signatures += 1
end

function basis_index(
    s :: Signature,
    signatures :: SignatureSet
) :: Int
    return signatures.basis_indices[s]
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

"""
Returns a list of all divisors of `s` in `signatures` that are distinct from
`s`.
"""
function enumerate_divisors(
    s :: Signature,
    signatures :: SignatureSet
) :: Vector{Signature}
    search_tree = signatures.signatures[s.index]
    return enumerate_reducers(s, Signature[], search_tree, skipbinomial=s)
end

#
# Definition of Signature Algorithms themselves.
#

const KoszulHeap{T} = BinaryHeap{Signature, ModuleMonomialOrdering{T}}
const SigLead = SignaturePolynomials.SigLead #I have no clue about why I need this...

mutable struct SignatureStats <: GBStats
    #These fields will be present in any GBStats type
    zero_reductions :: Int
    max_basis_size :: Int
    queued_pairs :: Int
    built_pairs :: Int
    reduced_pairs :: Int
    eliminated_by_truncation :: Int

    #These fields may be specific to the Signature algorithm
    eliminated_by_early_signature :: Int
    eliminated_by_late_signature :: Int
    eliminated_by_gcd :: Int
    eliminated_by_duplicate :: Int
    eliminated_by_early_koszul :: Int
    eliminated_by_late_koszul :: Int
    eliminated_by_base_divisors :: Int
    eliminated_by_low_base_divisors :: Int
    eliminated_by_high_base_divisors :: Int

    function SignatureStats()
        new(0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0)
    end
end

mutable struct SignatureAlgorithm{T, F <: Function} <: GBAlgorithm
    basis :: SigBasis{T}
    heap :: TriangleHeap{T, UInt32}
    koszul :: KoszulHeap{T}
    syzygies :: SignatureSet
    basis_signatures :: SignatureSet
    syzygy_pairs :: BitTriangle
    sigleads :: Vector{SigLead}
    comparator :: Comparator{SigLead, F}
    previous_signature :: Union{Signature, Nothing}
    should_truncate :: Bool
    stats :: SignatureStats
    preallocated :: Vector{Int}

    function SignatureAlgorithm(
        T :: Type,
        C :: Array{Float64, 2},
        A :: Array{Int, 2},
        b :: Vector{Int},
        mod_order :: Symbol,
        should_truncate :: Bool,
        minimization :: Bool
    )
        syzygies = SignatureSet()
        basis_signatures = SignatureSet()
        generators = SigPoly{T}[]
        sigleads = SigLead[]
        order = ModuleMonomialOrdering(C, A, b, mod_order, generators)
        module_order_lt = (s1, s2) -> SignaturePolynomials.signature_lt(s1, s2, order)
        basis = BinomialSet{SigPoly{T}, ModuleMonomialOrdering{T}}(generators, order, minimization)
        koszul = BinaryHeap{Signature}(order, [])
        syz_pairs = BitTriangle()
        comparator = Comparator{SigLead, typeof(module_order_lt)}(sigleads, module_order_lt)
        heap = TriangleHeap{T, UInt32}(basis, order, comparator)
        #Initialize all statistics collected by a SignatureAlgorithm
        stats = SignatureStats()
        new{T, typeof(module_order_lt)}(
            basis, heap, koszul, syzygies, basis_signatures, syz_pairs,
            sigleads, comparator, nothing, should_truncate, stats, Int[]
        )
    end
end

#
# Accessors and modifiers for SignatureAlgorithm
#

is_syzygy(i, j, algorithm :: SignatureAlgorithm) = algorithm.syzygy_pairs[i, j]
is_syzygy(pair :: SignaturePair, algorithm :: SignatureAlgorithm) =
    is_syzygy(pair.i, pair.j, algorithm)

function set_syzygy(i, j, algorithm :: SignatureAlgorithm)
    algorithm.syzygy_pairs[i, j] = true
end

set_syzygy(pair :: SignaturePair, algorithm :: SignatureAlgorithm) =
    set_syzygy(pair.i, pair.j, algorithm)

GBAlgorithms.use_implicit_representation(:: SignatureAlgorithm{T}) where {T} = is_implicit(T)

#
# GBAlgorithm interface implementation
#

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
    high_divisor :: Int,
    low_divisor :: Union{LowBaseDiv, Nothing},
    i :: Int,
    j :: Int
) :: Bool where {T <: GBElement}
    #GCD criterion
    if is_support_reducible(i, j, current_basis(algorithm))
        increment(algorithm, :eliminated_by_gcd)
        set_syzygy(i, j, algorithm)
        return true
    end
    if base_divisors_criterion(i, j, high_divisor, low_divisor, algorithm)
        increment(algorithm, :eliminated_by_base_divisors)
        set_syzygy(i, j, algorithm)
        return true
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
        increment(algorithm, :eliminated_by_early_signature)
        set_syzygy(pair, algorithm)
        return true
    end
    if koszul_criterion(pair, algorithm.koszul, algorithm.basis.order)
        increment(algorithm, :eliminated_by_early_koszul)
        set_syzygy(pair, algorithm)
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
        increment(algorithm, :eliminated_by_duplicate)
        set_syzygy(pair, algorithm)
        return true
    end
    algorithm.previous_signature = pair.signature
    ##Signature criterion
    #if signature_criterion(pair, algorithm.syzygies)
    #    data(algorithm)["eliminated_by_late_signature"] += 1
    #    set_syzygy(pair, algorithm)
    #    return true
    #end
    #if koszul_criterion(pair, algorithm.koszul, algorithm.basis.order)
    #    data(algorithm)["eliminated_by_late_koszul"] += 1
    #    set_syzygy(pair, algorithm)
    #    return true
    #end
    return false
end

"""
Adds `syzygy_element` to the set of known syzygies of `algorithm`.
"""
function GBAlgorithms.process_zero_reduction!(
    algorithm :: SignatureAlgorithm{T},
    syzygy_element :: SigPoly{T},
    pair :: SignaturePair
) where {T <: GBElement}
    add_signature!(algorithm.syzygies, signature(syzygy_element))
    set_syzygy(pair, algorithm)
    increment(algorithm, :zero_reductions)
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
    push!(current_basis(algorithm), copy(g))
    push!(algorithm.sigleads, g.siglead)
    FastComparator.update!(algorithm.comparator)
    add_signature!(algorithm.basis_signatures, g.signature, length(current_basis(algorithm)))
    add_row!(algorithm.syzygy_pairs)
    update_queue!(algorithm)
    #update_syzygies!(algorithm)
    #Add Koszul element to the queue if relevant
    if !isnothing(pair)
        k = koszul(GBElements.first(pair), GBElements.second(pair),
                   current_basis(algorithm), comparator=algorithm.comparator)
        if !isnothing(k)
            push!(algorithm.koszul, k)
        end
    end
    algorithm.stats.max_basis_size = max(algorithm.stats.max_basis_size,
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
    high_divisor = high_base_divisor(n, algorithm)
    low_divisor = low_base_divisor(n, algorithm)
    for i in 1:(n-1)
        if very_early_pair_elimination(algorithm, high_divisor, low_divisor, i, n)
            continue
        end
        sp = regular_spair(i, n, gb, comparator=algorithm.comparator)
        if !isnothing(sp) && !early_pair_elimination(algorithm, sp)
            push!(batch, sp)
            increment(algorithm, :queued_pairs)
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
        k = koszul(i, n, gb, comparator=algorithm.comparator)
        if !isnothing(k)
            push!(algorithm.koszul, k)
        end
    end
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
    C :: Array{Float64, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    # This assumes binary constraints
    # Must change the ordering due to new variables being added when C
    # was put into a standard form
    change_ordering!(current_basis(algorithm), C, A, b)
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    #Need to initialize the SignatureSets before adding new elements in case we
    #already start adding signature elements
    initialize!(algorithm.syzygies, num_gens)
    initialize!(algorithm.basis_signatures, num_gens)
    num_vars = size(A, 2)
    algorithm.preallocated = Vector{Int}(undef, num_vars)
    coef = zeros(Int, num_vars) #Coefficient of the signatures of these generators
    j = 1
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u, order(algorithm.basis),
                              check_truncation=truncate_basis(algorithm))
        if !isnothing(e)
            s = Signature(j, coef)
            GBAlgorithms.update!(algorithm, SigPoly{T}(e, s), nothing)
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
    if isequal(k_sig, spair.signature)
        return true
    end
    #TODO tiebreaking in the signature comparisons is wrong / isn't done
    #there may be bugs involving this elsewhere as well
    #if SignaturePolynomials.divides(k_sig, spair.signature)
    #    @show k_sig
    #    @show spair.signature
    #    @show Base.lt(order, k_sig, spair.signature)
    #    @show order.monomial_order * k_sig.monomial
    #    @show order.monomial_order * spair.signature.monomial
    #    @show order.monomial_order
    #    return true
    #end
    return false
end

"""
Applies the Base Divisors Criterion presented in Roune and Stillman (2012).
"""
function base_divisors_criterion(
    i :: Int,
    j :: Int,
    high_divisor :: Int,
    low_divisor :: Union{LowBaseDiv, Nothing},
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}
    if high_base_divisors_criterion(i, j, high_divisor, algorithm)
        increment(algorithm, :eliminated_by_high_base_divisors)
        return true
    end
    num_low_divisors = 2
    if low_base_divisors_criterion(i, j, low_divisor, algorithm)
        increment(algorithm, :eliminated_by_low_base_divisors)
        return true
    end
    return false
end

"""
Computes the index of the best high-ratio base divisor for the element of index
j.
"""
function high_base_divisor(
    j :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Int where {T <: GBElement}
    #Compute high-ratio base divisors
    gb = current_basis(algorithm)
    lt_b = leading_term(gb[j])
    base_divisors = enumerate_reducers(lt_b, gb, reduction_tree(gb), skipbinomial=gb[j])
    if isempty(base_divisors)
        return 0
    end
    #Find the one with lowest sig-lead ratio
    #TODO Get the indices of the base divisors from somewhere to avoid calling find_position
    comp = algorithm.comparator
    minimal_divisor = base_divisors[1]
    min_index = FastComparator.find_position(minimal_divisor.siglead, comp)
    #TODO I should pass this assert. Put it back in later.
    #@assert isequal(gb[min_index].siglead, minimal_divisor)
    for i in 2:length(base_divisors)
        div = base_divisors[i]
        div_index = FastComparator.find_position(div.siglead, comp)
        #TODO I should pass this assert. Put it back in later.
        #@assert isequal(gb[div_index].siglead, div)
        if FastComparator.compare(comp, div_index, min_index) == :lt
            minimal_divisor = div
            min_index = div_index
        end
    end
    return min_index
end

function low_base_divisor(
    j :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Union{LowBaseDiv, Nothing} where {T <: GBElement}
    gb = current_basis(algorithm)
    comp = algorithm.comparator
    sig_b = gb[j].signature
    base_sigs = enumerate_divisors(sig_b, algorithm.basis_signatures)
    if isempty(base_sigs)
        return nothing
    end
    max_index = basis_index(base_sigs[1], algorithm.basis_signatures)
    max_div = gb[max_index]
    for i in 2:length(base_sigs)
        div_index = basis_index(base_sigs[i], algorithm.basis_signatures)
        div = gb[div_index]
        if FastComparator.compare(comp, div_index, max_index) == :gt
            max_index = div_index
            max_div = div
        end
    end
    return LowBaseDiv(max_index, max_div, gb[j])
end

"""
The Base divisor criterion for high-ratio elements. An element a is high-ratio
if lt(a) | lt(b), where b is the new element being added.
"""
function high_base_divisors_criterion(
    i :: Int,
    j :: Int,
    divisor_index :: Int,
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}
    if divisor_index == 0 #There's no high-ratio divisor
        return false
    end
    #There's a high-ratio divisor, attempt to use it
    if i != divisor_index && is_syzygy(i, divisor_index, algorithm)
        if FastComparator.compare(algorithm.comparator, i, divisor_index) == :gt &&
            FastComparator.compare(algorithm.comparator, i, j) == :gt
            return true
        end
    end
    return false
end

"""
The Base divisor criterion for low-ratio elements. An element a is low-ratio if
signature(a) | signature(b), where b is the new element being added.
"""
function low_base_divisors_criterion(
    i :: Int,
    j :: Int,
    low_divisor :: Union{LowBaseDiv, Nothing},
    algorithm :: SignatureAlgorithm{T}
) :: Bool where {T <: GBElement}
    if isnothing(low_divisor) #There's no low-ratio divisor
        return false
    end
    #There's a low-ratio divisor, attempt to use it
    comp = algorithm.comparator
    if i != low_divisor.index && is_syzygy(i, low_divisor.index, algorithm)
        if FastComparator.compare(comp, i, low_divisor.index) == :lt &&
            FastComparator.compare(comp, i, j) == :lt
            lt_c = leading_term(current_basis(algorithm)[i])
            if all(i -> lt_c[i] <= low_divisor.monomial[i], 1:length(lt_c))
                return true
            end
        end
    end
    return false
end

end
