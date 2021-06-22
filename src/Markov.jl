"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis

TODO implement a Project-and-Lift state structure. This will facilitate the
implementation of group relaxations and early termination, as well as
iteratively computing parts of the GB and continuing later
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
TODO I probably could make projections more efficient. But at least
this way they are very simple to implement...
"""
function project_and_lift(
    instance :: IPInstance;
    truncate :: Bool = true
) :: Vector{Vector{Int}}
    #Compute lattice basis, TODO refactor
    m, n = size(instance.A)
    Mmn = MatrixSpace(ZZ, m, n)
    rank, basis = kernel(Mmn(instance.A))
    basis = basis' #Put the basis in row form
    hnf_basis = hnf(basis)
    normalize_hnf!(hnf_basis)
    @show hnf_basis
    #Now, the first `rank` columns of hnf_basis should be LI
    #and thus a Markov basis of the corresponding projection
    sigma = Vector(rank+1:n)
    nonnegative = nonnegative_vars(instance)
    for s in sigma
        nonnegative[s] = false
    end
    #This is a Markov basis of the projected problem
    markov = Vector{Int}[]
    for row in eachrow(Array(hnf_basis))
        v = Vector{Int}(row)
        push!(markov, v)
    end
    projection = IPInstance(instance.A, instance.b, instance.C, instance.u,
                            nonnegative)
    while !isempty(sigma)
        i = random_lifting(sigma) #Pick some variable to lift
        perm_i = projection.permutation[i] #This is the index of i in projection
        if is_bounded(perm_i, projection)
            #Hard case, compute GB
        else
            #Easy case
        end
        #Finished lifting i, remove it from sigma
        i_index = findfirst(isequal(i), sigma)
        deleteat!(sigma, i_index)
        nonnegative[i] = true
        projection = IPInstance(instance.A, instance.b, instance.C, instance.u,
                                nonnegative)
    end
    return markov
end

"""
Computes a Markov basis with the simplified algorithm. Assumes all data in
instance.A and instance.b is non-negative.
"""
function simple_markov(
    instance :: IPInstance
) :: Vector{Vector{Int}}
    @assert is_nonnegative(instance) #Check whether the hypotheses hold
    n = size(instance.A, 2) - size(instance.A, 1)
    m = size(instance.A, 1) - n
    basis = Vector{Int}[]
    for i in 1:n
        v = zeros(Int, n)
        v[i] = 1
        r = zeros(Int, n)
        r[i] = -1
        s = -copy(instance.A[1:m, i])
        g = vcat(v, s, r)
        push!(basis, g)
    end
    return basis
end

"""
Computes a Markov basis of `instance`.
"""
function markov_basis(
    instance :: IPInstance;
    algorithm :: Symbol = :Any
) :: Vector{Vector{Int}}
    if algorithm == :Any
        if is_nonnegative(instance)
            #The Simple algorithm may be used, so use it.
            algorithm = :Simple
        else
            #The Simple algorithm can't be used, so fall back to project-and-lift
            algorithm = :ProjectAndLift
        end
    end
    if algorithm == :Simple
        basis = simple_markov(instance)
    else
        basis = project_and_lift(instance)
    end
    return basis
end

"""
Computes a Markov basis for the IP instance

min C * x
s.t. A * x = b
0 <= x <= u
"""
function markov_basis(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    algorithm :: Symbol = :Any
) :: Vector{Vector{Int}}
    instance = IPInstance(A, b, C, u)
    return markov_basis(instance, algorithm=algorithm)
end

end
