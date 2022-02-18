module FGLM

using IPGBs.BinomialSets
using IPGBs.GBElements
using IPGBs.Orders
using IPGBs.MonomialHeaps

"""
    update_monomial_heap!(h :: MonomialHeap{T}, m :: Vector{Int}) where {T <: GBOrder}

Add all monomials of the form x_i * `m` for all variables x_i to the heap `h` and count repetitions.
"""
function update_monomial_heap!(
    h :: MonomialHeap{S},
    m :: Vector{Int},
    nf :: T
) where {T, S <: GBOrder}
    n = length(m)
    for i in 1:n
        new_monom = copy(m)
        new_monom[i] += 1
        push!(h, new_monom, nf, i) #Push will deal with the repetitions
    end
end

"""
    find_linear_dependency(nf :: Vector{Int}, std_basis :: Dict{Vector{Int}, Vector{Int}}

Return an element of `std_basis` with normal form `nf`, if such an element
exists, otherwise return [].

TODO It is likely possible to do this more efficient using some special kind
of tree. I'll leave this for later, in case it's needed.
"""
function find_linear_dependency(
    nf :: Vector{Int},
    std_basis :: Dict{Vector{Int}, Vector{Int}}
) :: Vector{Int}
    if haskey(std_basis, nf)
        return std_basis[nf]
    end
    return []
end

"""
    is_below_staircase(m :: Vector{Int}, gb :: Vector{T}) where {T}

Checks whether `m` is below the staircase given by `gb`.
Naïve implementation. It is possible to do this more efficiently by counting
the number of times FGLM generated `m`.
"""
function is_below_staircase(
    m :: WeightedMonomial,
    gb :: Vector{T}
) :: Bool where {T}
    for basis_monom in gb
        if all(basis_monom[i] <= m[i] for i in 1:length(m))
            return false
        end
    end
    return true
end

"""
    is_below_staircase_fast(m :: WeightedMonomial)

Check whether `m` is below a GB staircase indirectly by counting how many times
it was inserted in the priority queue.

This idea comes straight from the original FGLM paper Faugère et al (1994).
"""
function is_below_staircase_fast(
    m :: WeightedMonomial
) :: Bool
    return m.count < length(m.monomial)
end

"""
    direct_normal_form(m :: WeightedMonomial, gb :: BinomialSet{T, S}, target_order :: S) where {T, S <: GBOrder}

Compute the normal form of `m` with respect to `gb` by applying the usual reduction process directly.
"""
function direct_normal_form(
    m :: WeightedMonomial,
    gb :: BinomialSet{T, S},
    target_order :: S
) :: T where {T, S <: GBOrder}
    m' = to_gbelement(copy(m), target_order, T)
    nf = BinomialSets.reduce!(m', gb)
    return nf
end

"""
    fast_normal_form(m :: WeightedMonomial, gb :: BinomialSet{T, S}, target_order :: S) where {T, S <: GBOrder}

Compute the normal form of `m` with respect to `gb` by computing NormalForm(x_i * NormalForm(m')) where m = x_i * m'.

This is the method suggested in the original FGLM paper Faugère et al (1994).
"""
function fast_normal_form(
    m :: WeightedMonomial,
    gb :: BinomialSet{T, S},
    target_order :: S
) :: T where {T, S <: GBOrder}
    #Two cases: either m is 1, so there is no previously known normal form,
    #or m isn't 1, and thus we have a divisor and var index in m.
    if m.var_index == 0
        #Case 1: no previous normal form
        return direct_normal_form(m, gb, target_order)
    end
    #Case 2: use the divisor's normal form to speed up the computation
    prev_nf = copy(m.divisor_nf)
    prev_nf[m.var_index] += 1
    nf = BinomialSets.reduce!(prev_nf, gb)
    return nf
end

"""
    fglm(gb1 :: BinomialSet{T, S}, target_order :: S) where {T, S <: GBOrder}

Convert a Gröbner basis `gb1` to `target_order` using the FGLM algorithm.

FGLM only works when the ideal generated by `gb1` is zero-dimensional. This
is assumed in the implementation.
"""
function fglm(
    gb1 :: BinomialSet{T, S},
    target_order :: S
) :: BinomialSet{T, S} where {T, S <: GBOrder}
    if isempty(gb1)
        return BinomialSet(T[], target_order)
    end
    n = length(gb1[1])
    one = zeros(Int, n) #The monomial 1
    next_monomials = MonomialHeap(target_order, [ one ])
    #The std_basis dict maps normal forms w.r.t the source order to standard
    #monomials in the target order
    std_basis = Dict{Vector{Int}, Vector{Int}}()
    gb2 = T[]
    while !isempty(next_monomials)
        m = pop!(next_monomials)
        if is_below_staircase_fast(m)
            nf = fast_normal_form(m, gb1, target_order)
            ld = find_linear_dependency(nf, std_basis)
            if isempty(ld)
                #Linearly independent case, new std_basis monomial
                std_basis[nf] = m.monomial
                update_monomial_heap!(next_monomials, m.monomial, nf)
            else
                #Linearly dependent case, new gb2 binomial
                new_binomial = m.monomial - ld
                new_elem = to_gbelement(new_binomial, target_order, T)
                push!(gb2, new_elem)
            end
        end
    end
    return BinomialSet(gb2, target_order)
end

end
