"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis
"""
module Markov

export markov_basis

using AbstractAlgebra

using IPGBs.IPInstances

"""
Changes `H` to an upper row HNF matrix satisfying the following property:
- all entries above a pivot are non-positive and of smaller magnitude than the pivot

Assumes `H` is already in HNF as defined by AbstractAlgebra.jl, that is, it is in
upper row HNF satisfying:
- all entries above a pivot are non-negative and of smaller magnitude than the pivot
"""
function normalize_hnf!(
    H :: Generic.MatSpaceElem{T}
) where {T}
    m, n = size(H)
    for i in 1:m
        #Find the pivot in row i
        j = i + 1
        while j <= n && H[i, j] == 0
            j += 1
        end
        if j > n #We reached a row of zeros, we are done.
            break
        end
        #Update rows above the pivot
        for k in 1:(i-1)
            if H[k, j] > 0 #only change positive entries
                H[k, :] -= H[i, :]
            end
        end
    end
end

"""
Verifies whether `H` is in normalized HNF form as defined in `normalize_hnf!`.
"""
function is_normalized_hnf(
    H :: Generic.MatSpaceElem{T}
) :: Bool where {T}
    m, n = size(H)
    for i in 1:m
        j = i + 1
        while j <= n && H[i, j] == 0
            j += 1
        end
        if j > n
            break
        end
        for k in 1:(i-1)
            if H[k, j] > 0 || (H[k, j] < 0 && abs(H[k, j]) >= H[i, j])
                return false
            end
        end
    end
    return true
end

"""
Given a saturated lattice in the form L = {x | Ax = 0} find a lattice basis
for it using the AbstractAlgebra library.
"""
function lattice_basis(
    A :: Array{Int, 2}
) :: Vector{Vector{BigInt}}
    m, n = size(A)
    space = MatrixSpace(ZZ, m, n)
    rank, basis = kernel(space(A))
    julia_form = Array(basis)
    return [ julia_form[:, j] for j in 1:size(julia_form, 2)]
end

random_lifting(sigma) = rand(sigma)
minimum_lifting(sigma) = minimum(sigma)

"""
TODO understand what exactly will be the inputs of this function
Is it the original problem formulation or the lattice basis or...?
"""
function project_and_lift(
    instance :: IPInstance
) :: Vector{Vector{BigInt}}
    #Compute lattice basis, TODO refactor
    m, n = size(instance.A)
    Mmn = MatrixSpace(ZZ, m, n)
    rank, basis = kernel(Mmn(instance.A))
    basis = basis' #Put the basis in row form
    hnf_basis = hnf(basis)
    normalize_hnf!(hnf_basis)
    #Now, the first `rank` columns of hnf_basis should be LI
    #and thus a Markov basis of the corresponding projection
    #TODO use that to compute sigma below
    sigma = []
    while !isempty(sigma)
        i = random_lifting(sigma) #Pick some variable to lift
        if unbounded(i)
            #Easy case
        else
            #Hard case, compute GB
        end
    end
end

"""
Computes a Markov basis of `instance`.
"""
function markov_basis(
    instance :: IPInstance,
    algorithm = :Any
)
    #Do preprocessing, select algorithm, etc...
    project_and_lift(instance)
end

end
