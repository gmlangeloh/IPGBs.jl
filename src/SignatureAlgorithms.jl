"""
An implementation of a Signature-based algorithm for Gröbner bases of toric
ideals.
"""
module SignatureAlgorithms

using DataStructures

using IPGBs.GBElements
using IPGBs.SupportTrees
using IPGBs.SignaturePolynomials

"""
Returns the generators used for a truncated GB of a toric ideal given by an IP,
assuming all data is non-negative.
"""
function truncated_generators(
    A :: Array{Int, 2},
    C :: Array{Int, 2}
) :: Vector{SigPoly}
    generators = Vector{SigPoly}()
    n = size(A, 2)
    coef = zeros(Int, n) #Coefficient of the signatures of these generators
    for i in 1:n
        e = lattice_generator(i, A, C)
        s = Signature(i, coef)
        push!(generators, SigPoly(e, s))
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
    generators :: SigBasis,
    module_ordering :: ModuleMonomialOrdering
) #TODO type this, BinaryHeap...
    heap = BinaryHeap{SPair}(module_ordering, [])
    for i in 1:length(generators)
        for j in 1:(i-1)
            sp = regular_spair(i, j, generators)
            if sp != nothing
                push!(heap, sp)
            end
        end
    end
    return heap
end

function update_queue!(
    spairs,
    gb :: SigBasis
)
    n = length(gb)
    for i in 1:(n-1)
        sp = regular_spair(i, n, gb)
        if sp != nothing
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
            #println("Eliminating by sig criterion")
            #@show spair.i spair.j spair.signature syz
            return true
        end
    end
    return false
end

"""
Returns true iff this spair is eliminated by some criterion.
"""
function post_criteria(
    spair :: SPair,
    syzygies :: Vector{Signature}
) :: Bool
    #TODO add other criteria as necessary
    return signature_criterion(spair, syzygies)
end

function signature_algorithm(
    generators :: Vector{SigPoly},
    order :: Array{Int, 2},
    module_ordering :: ModuleMonomialOrdering,
    A :: Array{Int, 2},
    b :: Vector{Int},
    u :: Vector{Int}
) :: Vector{Vector{Int}}
    println("STARTING")
    reducer = support_tree(generators, fullfilter=true)
    gb = SigBasis(copy(generators), module_ordering, reducer)
    spairs = make_priority_queue(gb, module_ordering)
    syzygies = Vector{Signature}() #TODO maybe this should be some other
    #structure instead of a vector... Also, maybe I should generate the
    #Koszul syzygies initially.
    reduction_count = 0
    previous_sig = nothing
    while !isempty(spairs)
        sp = pop!(spairs)
        #Process each signature only once
        if previous_sig != nothing && isequal(previous_sig, sp.signature)
            continue
        end
        previous_sig = sp.signature
        if post_criteria(sp, syzygies)
            #eliminated = build_spair(sp, gb)
            #@show eliminated
            continue
        end
        p = build_spair(sp, gb)
        if !isfeasible(p.polynomial, A, b, u)
            continue
        end
        println("before")
        @show length(gb)
        @show sp.i sp.j
        @show p p.signature
        SupportTrees.reduce!(p, gb, gb.reduction_tree)
        reduction_count += 1
        println("after")
        @show p p.signature
        if SignaturePolynomials.iszero(p)
            println("new syzygy")
            @show p.signature
            push!(syzygies, p.signature)
        else
            push!(gb, p)
            #println("new basis")
            update_queue!(spairs, gb)
        end
    end
    for i in 1:length(gb)
        g = gb[i]
        println(i, " ", g.polynomial, " ", g.signature)
    end
    #Return the representation of each element as an integer vector
    @show reduction_count
    return projection.(gb)
end

function siggb(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    module_order :: ModuleMonomialOrder = SignaturePolynomials.ltpot
) :: Vector{Vector{Int}}
    order = make_monomial_order(C)
    generators = truncated_generators(A, C)
    #We can check the comparisons between signatures
    #for i in 1:length(generators)
    #    for j in 1:(i-1)
    #        si = generators[i].signature
    #        sj = generators[j].signature
    #        r = SignaturePolynomials.ltpot_lt(si, sj, order, generators)
    #        println(i, " < ", j, " = ", r)
    #    end
    #end
    sig_ordering = ModuleMonomialOrdering(C, module_order, generators)
    return signature_algorithm(
        generators, order, sig_ordering, A, b, u
    )
end

end
