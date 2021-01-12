"""
An implementation of a Signature-based algorithm for Gröbner bases of toric
ideals.
"""
module SignatureAlgorithms
export siggb

using DataStructures

using IPGBs.Buchberger
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.Binomials
using IPGBs.GradedBinomials
using IPGBs.SupportTrees
using IPGBs.SignaturePolynomials

"""
Returns the generators used for a truncated GB of a toric ideal given by an IP,
assuming all data is non-negative.
"""
function truncated_generators(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int},
    T :: DataType,
    lattice_generator :: Function
) :: Vector{SigPoly{T}}
    generators = Vector{SigPoly{T}}()
    # This assumes binary constraints
    if T == Binomial
        num_gens = size(A, 2) - size(A, 1)
    else
        num_gens = size(A, 2)
    end
    num_vars = size(A, 2)
    coef = zeros(Int, num_vars) #Coefficient of the signatures of these generators
    j = 1
    for i in 1:num_gens
        e = lattice_generator(i, A, b, C, u)
        if !isnothing(e)
            s = Signature(j, coef)
            push!(generators, SigPoly{T}(e, s))
            j += 1
        end
    end
    return generators
end

"""
Adds a tiebreaking order to the given monomial order `C`.
"""
function make_monomial_order(
    C :: Array{Int, 2};
    tiebreaker :: String = "grevlex"
) :: Array{Int, 2}
    n = size(C, 2)
    if tiebreaker == "grevlex"
        tie_matrix = [ i <= (n+1) - j ? 1 : 0 for i=1:n, j=1:n ]
    else #Use lex tiebreaker
        tie_matrix = [ i == j ? 1 : 0 for i=1:n, j=1:n ]
    end
    return vcat(C, tie_matrix)
end

"""
Create priority queue and add all regular S-pairs built from generators to it.
"""
function make_priority_queue(
    generators :: SigBasis{T},
    module_ordering :: ModuleMonomialOrdering
) where {T <: GBElement}
    #TODO type this, BinaryHeap...
    heap = BinaryHeap{SPair}(module_ordering, [])
    for i in 1:length(generators)
        for j in 1:(i-1)
            sp = regular_spair(i, j, generators)
            if !isnothing(sp)
                push!(heap, sp)
            end
        end
    end
    return heap
end

function update_queue!(
    spairs, #TODO type this thing. It is supposed to be the heap
    gb :: SigBasis{T}
) where {T <: GBElement}
    n = length(gb)
    for i in 1:(n-1)
        sp = regular_spair(i, n, gb)
        if !isnothing(sp)
            push!(spairs, sp)
        end
    end
end

"""
Returns true iff spair is eliminated by the signature criterion, i.e. there
is a known syzygy that divides its signature.
"""
function signature_criterion(
    spair :: SPair,
    syzygies :: Vector{Signature}
) :: Bool
    for syz in syzygies
        if divides(syz, spair.signature)
            return true
        end
    end
    return false
end

"""
Returns true iff this spair is eliminated by some criterion.
"""
function late_criteria(
    spair :: SPair,
    gb :: SigBasis{T},
    syzygies :: Vector{Signature},
    minimization :: Bool
) :: Bool where {T <: GBElement}
    if signature_criterion(spair, syzygies)
        return true
    end
    if Buchberger.is_support_reducible( #GCD criterion
        spair.i, spair.j, gb.positive_supports, gb.negative_supports, minimization
    )
        return true
    end
    return false
end

"""
Creates a list of all Koszul syzygies of gb.
"""
function initial_syzygies(
    gb :: SigBasis{T}
) :: Vector{Signature} where {T <: GBElement}
    syzygies = Vector{Signature}()
    for i in 1:length(gb)
        for j in 1:(i-1)
            push!(syzygies, koszul(i, j, gb))
        end
    end
    return syzygies
end

"""
Adds the Koszul syzygies corresponding to the newest element of gb to the
syzygy list.
"""
function update_syzygies!(
    syzygies :: Vector{Signature},
    gb :: SigBasis{T}
) where {T <: GBElement}
    n = length(gb)
    for i in 1:(n-1)
        push!(syzygies, koszul(i, n, gb))
    end
end

function signature_algorithm(
    generators :: Vector{SigPoly{T}},
    order :: Array{Int, 2},
    module_ordering :: ModuleMonomialOrdering,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int},
    structure :: DataType,
    minimization :: Bool
) :: Vector{Vector{Int}} where {T <: GBElement}
    reducer = support_tree(generators, fullfilter=(structure == GradedBinomial))
    gb = SigBasis(copy(generators), module_ordering, reducer)
    spairs = make_priority_queue(gb, module_ordering)
    syzygies = initial_syzygies(gb)
    reduction_count = 0
    previous_sig = nothing
    while !isempty(spairs)
        sp = pop!(spairs)
        #Process each signature only once
        if !isnothing(previous_sig) && isequal(previous_sig, sp.signature)
            continue
        end
        previous_sig = sp.signature
        if late_criteria(sp, gb, syzygies, minimization)
            continue
        end
        p = build_spair(sp, gb)
        if !isfeasible(p, A, b, u)
            continue
        end
        reduced_to_zero = SupportTrees.reduce!(p, gb, gb.reduction_tree)
        reduction_count += 1
        if reduced_to_zero
            push!(syzygies, p.signature)
        else
            push!(gb, p)
            update_queue!(spairs, gb)
            update_syzygies!(syzygies, gb)
        end
    end
    #Return the representation of each element as an integer vector
    if structure == Binomial
        output_basis = projection.(gb)
    else
        output_basis = [ -GradedBinomials.fullform(g.polynomial) for g in gb ]
    end
    @show reduction_count
    @show length(syzygies)
    return output_basis
end

function initial_gb(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    T :: DataType = Binomial,
) :: Vector{SigPoly{T}}
    if T == Binomial
        lattice_generator = lattice_generator_binomial
    else
        lattice_generator = lattice_generator_graded
    end
    generators = truncated_generators(
        A, b, C, u, T, lattice_generator
    )
    return generators
end

function siggb(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    module_order :: ModuleMonomialOrder = SignaturePolynomials.ltpot,
    structure :: DataType = Binomial
) :: Vector{Vector{Int}}
    order = make_monomial_order(C)
    minimization = structure == Binomial
    A, b, C, u = GBTools.normalize(
        A, b, C, u, apply_normalization=minimization
    )
    generators = initial_gb(A, b, C, u, T = structure)
    sig_ordering = ModuleMonomialOrdering(C, module_order, generators)
    return signature_algorithm(
        generators, order, sig_ordering, A, b, u, structure, minimization
    )
end

end
