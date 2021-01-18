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

using IPGBs.GBAlgorithms

#
# Definition of Signature Algorithms themselves.
#

struct SignatureAlgorithm{T} <: GBAlgorithm
    basis :: SigBasis{T}
    heap #TODO type this!!!
    syzygies :: Vector{Signature}
    previous_signature :: Union{Signature, Nothing}
    stats :: GBStats

    function SignatureAlgorithm(T :: Type, C :: Array{Int, 2})
        syzygies = Signature[]
        generators = SigPoly{T}[]
        #For now, fix the Schreyer order as module monomial order
        order = ModuleMonomialOrdering(C, SignaturePolynomials.ltpot, generators)
        basis = BinomialSet(generators, order)
        heap = BinaryHeap{SignaturePair}(order, [])
        #TODO create additional stats fields for Signatures
        #This should include counting the number of pairs eliminated by each
        #criterion
        new{T}(basis, heap, syzygies, nothing, GBStats())
    end
end

stats(algorithm :: SignatureAlgorithm) = algorithm.stats.stats
current_basis(algorithm :: SignatureAlgorithm) = algorithm.basis

function next_pair!(
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
    #TODO implement this!
    return false
end

function late_pair_elimination(
    algorithm :: SignatureAlgorithm{T},
    pair :: SignaturePair
) :: Bool where {T <: GBElement}
    #Skip repeated signatures
    if !isnothing(algorithm.previous_signature) &&
        isequal(algorithm.previous_signature, pair.signature)
        return true
    end
    algorithm.previous_signature = pair.signature
    #Signature criterion
    if signature_criterion(pair, algorithm.syzygies)
        return true
    end
    #GCD criterion
    if is_support_reducible(first(pair), second(pair), current_basis(algorithm))
        return true
    end
    return false
end

function process_zero_reduction(
    algorithm :: GBAlgorithm{T},
    syzygy_element :: SigPoly{T}
) where {T <: GBElement}
    push!(algorithm.syzygies, signature(syzygy_element))
    stats(algorithm)["zero_reductions"] += 1
end

function update!(
    algorithm :: GBAlgorithm,
    g :: SigPoly{T}
) where {T <: GBElement}
    push!(current_basis(algorithm), g)
    update_queue!(algorithm)
    update_syzygies!(algorithm)
    stats(algorithm)["max_basis_size"] = max(stats(algorithm)["max_basis_size"],
                                             length(current_basis(algorithm)))
end

"""
Adds new SignaturePairs to algorithm.heap. Applies early elimination of S-pairs
when possible.
"""
function update_queue!(algorithm :: SignatureAlgorithm)
    gb = current_basis(algorithm)
    n = length(gb)
    for i in 1:(n-1)
        sp = regular_spair(i, n, gb)
        if !isnothing(sp) && !early_pair_elimination(algorithm, sp)
            push!(algorithm.heap, sp)
            stats(algorithm)["queued_pairs"] += 1
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
        push!(algorithm.syzygies, koszul(i, n, gb))
    end
end

#TODO almost all of this logic is still duplicated in Buchberger.jl. Can I do
#this better?
function initialize!(
    algorithm :: SignatureAlgorithm{T},
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int}
) where {T <: GBElement}
    # This assumes binary constraints
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
        lattice_generator = lattice_generator_binomial
    else
        num_gens = size(A, 2)
        lattice_generator = lattice_generator_graded
    end
    num_vars = size(A, 2)
    coef = zeros(Int, num_vars) #Coefficient of the signatures of these generators
    j = 1
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u)
        if !isnothing(e)
            s = Signature(j, coef)
            update!(algorithm, SigPoly{T}(e, s))
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
    syzygies :: Vector{Signature}
) :: Bool
    for syz in syzygies
        if divides(syz, spair.signature)
            return true
        end
    end
    return false
end

#TODO the remainder of this module has to be refactored

#"""
#Returns the generators used for a truncated GB of a toric ideal given by an IP,
#assuming all data is non-negative.
#
#TODO slightly change the arguments here so that it coincides with the other
#generators. Then define a initial_generators for any algorithm, and call functions
#like this one depending on the type that is inside
#"""
#function truncated_generators(
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    C :: Array{Int, 2},
#    u :: Vector{Int},
#    T :: DataType,
#    lattice_generator :: Function
#) :: Vector{SigPoly{T}}
#    generators = Vector{SigPoly{T}}()
#    # This assumes binary constraints
#    if T == Binomial
#        num_gens = size(A, 2) - size(A, 1)
#    else
#        num_gens = size(A, 2)
#    end
#    num_vars = size(A, 2)
#    coef = zeros(Int, num_vars) #Coefficient of the signatures of these generators
#    j = 1
#    for i in 1:num_gens
#        e = lattice_generator(i, A, b, C, u)
#        if !isnothing(e)
#            s = Signature(j, coef)
#            push!(generators, SigPoly{T}(e, s))
#            j += 1
#        end
#    end
#    return generators
#end

#"""
#Adds a tiebreaking order to the given monomial order `C`.
#
#TODO this thing should definitely be somewhere else
#"""
#function make_monomial_order(
#    C :: Array{Int, 2};
#    tiebreaker :: String = "grevlex"
#) :: Array{Int, 2}
#    n = size(C, 2)
#    if tiebreaker == "grevlex"
#        tie_matrix = [ i <= (n+1) - j ? 1 : 0 for i=1:n, j=1:n ]
#    else #Use lex tiebreaker
#        tie_matrix = [ i == j ? 1 : 0 for i=1:n, j=1:n ]
#    end
#    return vcat(C, tie_matrix)
#end

#"""
#Create priority queue and add all regular S-pairs built from generators to it.
#"""
#function make_priority_queue(
#    generators :: SigBasis{T}
#) where {T <: GBElement}
#    #TODO type this, BinaryHeap...
#    heap = BinaryHeap{SPair}(order(generators), [])
#    for i in 1:length(generators)
#        for j in 1:(i-1)
#            sp = regular_spair(i, j, generators)
#            if !isnothing(sp)
#                push!(heap, sp)
#            end
#        end
#    end
#    return heap
#end

#"""
#Creates a list of all Koszul syzygies of gb.
#
#TODO can be done in the initialization function
#"""
#function initial_syzygies(
#    gb :: SigBasis{T}
#) :: Vector{Signature} where {T <: GBElement}
#    syzygies = Vector{Signature}()
#    for i in 1:length(gb)
#        for j in 1:(i-1)
#            push!(syzygies, koszul(i, j, gb))
#        end
#    end
#    return syzygies
#end
#
#function signature_algorithm(
#    generators :: Vector{SigPoly{T}},
#    order :: Array{Int, 2},
#    module_ordering :: ModuleMonomialOrdering,
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    u :: Vector{Int},
#    structure :: DataType,
#    minimization :: Bool
#) :: Vector{Vector{Int}} where {T <: GBElement}
#    gb :: SigBasis{T} = BinomialSet(generators, module_ordering)
#    spairs = make_priority_queue(gb)
#    syzygies = initial_syzygies(gb)
#    reduction_count = 0
#    zero_reductions = 0
#    previous_sig = nothing
#    while !isempty(spairs)
#        sp = pop!(spairs)
#        #Process each signature only once
#        if !isnothing(previous_sig) && isequal(previous_sig, sp.signature)
#            continue
#        end
#        previous_sig = sp.signature
#        if late_criteria(sp, gb, syzygies)
#            continue
#        end
#        p = build_spair(sp, gb)
#        if !isfeasible(p, A, b, u)
#            continue
#        end
#        reduced_to_zero = BinomialSets.reduce!(p, gb)
#        reduction_count += 1
#        if reduced_to_zero
#            push!(syzygies, p.signature)
#            zero_reductions += 1
#        else
#            push!(gb, p)
#            update_queue!(spairs, gb)
#            update_syzygies!(syzygies, gb)
#        end
#    end
#    #Return the representation of each element as an integer vector
#    output_basis = BinomialSets.fourti2_form(gb)
#    @show reduction_count zero_reductions
#    @show length(syzygies)
#    return output_basis
#end
#
#function initial_gb(
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    C :: Array{Int, 2},
#    u :: Vector{Int};
#    T :: DataType = Binomial,
#) :: Vector{SigPoly{T}}
#    if T == Binomial
#        lattice_generator = lattice_generator_binomial
#    else
#        lattice_generator = lattice_generator_graded
#    end
#    generators = truncated_generators(
#        A, b, C, u, T, lattice_generator
#    )
#    return generators
#end
#
#function siggb(
#    A :: Array{Int, 2},
#    b :: Vector{Int},
#    C :: Array{Int, 2},
#    u :: Vector{Int};
#    module_order :: ModuleMonomialOrder = SignaturePolynomials.ltpot,
#    structure :: DataType = Binomial
#) :: Vector{Vector{Int}}
#    order = make_monomial_order(C)
#    minimization = structure == Binomial
#    A, b, C, u = GBTools.normalize(
#        A, b, C, u, apply_normalization=minimization
#    )
#    generators = initial_gb(A, b, C, u, T = structure)
#    sig_ordering = ModuleMonomialOrdering(C, module_order, generators)
#    return signature_algorithm(
#        generators, order, sig_ordering, A, b, u, structure, minimization
#    )
#end

end
