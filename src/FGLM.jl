module FGLM

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Orders

"""
Auxiliary struct for FGLM. Represents a monomial in the standard basis of the
target order in the algorithm. For efficiency, these monomials are stored
alongside their normal forms, which are known.
"""
struct StdBasisMonomial
    monomial :: Vector{Int}
    normal_form :: Vector{Int}
end

"""
Returns 0 if there's no element in `std_basis` with normal form `nf`.
Otherwise, returns the index of an element with this normal form.

TODO this is terribly inefficient. This search can be done using a tree or
in some other way.
"""
function find_linear_dependency(
    nf :: Vector{Int},
    std_basis :: Vector{StdBasisMonomial}
) :: Int
    for i in 1:length(std_basis)
        if std_basis[i].normal_form == nf
            return i
        end
    end
    return 0
end

function update_next_monomials!(
    next_monomials :: Vector{Vector{Int}},
    m :: Vector{Int},
    order :: S
) where {S <: GBOrder}
    #Add new monomials to next_monomials, keeping it ordered and avoiding
    #repetitions
end

function fglm(
    gb1 :: BinomialSet{T, S},
    target_order :: S
) :: BinomialSet{T, S} where {T, S <: GBOrder}
    if isempty(gb1)
        return BinomialSet(T[], target_order)
    end
    n = length(gb1[1])
    one = zeros(Int, n) #The monomial 1
    next_monomials = [ one ]
    std_basis = StdBasisMonomial[]
    gb2 = T[]
    while !isempty(next_monomials)
        m = popfirst!(next_monomials)
        if m #TODO Isn't a multiple of some previous lt, implement test
            normal_form = BinomialSets.reduce!(m, gb1) #TODO the typing in reduce probably doesn't allow this, I'll have to fix it later
            ld = find_linear_dependency(normal_form, std_basis)
            if ld == 0 #Linearly independent case, new std_basis monomial
                push!(std_basis, StdBasisMonomial(m, normal_form))
                update_next_monomials!(next_monomials, m, target_order)
            else #Linearly dependent case, new gb2 binomial
                new_binomial = m - std_basis[ld].monomial
                new_elem = to_gbelement(T, new_binomial, target_order)
                push!(gb2, new_elem)
            end
        end
    end
    return BinomialSet(gb2, target_order)
end

end
