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

using IPGBs.Binomials
using IPGBs.Buchberger
using IPGBs.GBAlgorithms
using IPGBs.GBElements
using IPGBs.IPInstances

"""
    normalize_hnf!(H :: Generic.MatSpaceElem{T})

Change `H` to an upper row HNF matrix satisfying the following property:
all entries above a pivot are non-positive and of smaller magnitude than the pivot

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
    is_normalized_hnf(H :: Generic.MatSpaceElem{T})

Return true iff `H` is in normalized HNF form as defined in `normalize_hnf!`.
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
    lattice_basis(A :: Matrix{Int}) :: Vector{Vector{BigInt}}

Find a lattice basis for the saturated lattice in the form L = {x | Ax = 0}
using the AbstractAlgebra library.

TODO Do I need this? Why? Review project-and-lift.
"""
function lattice_basis(
    A :: Matrix{Int}
) :: Vector{Vector{BigInt}}
    m, n = size(A)
    space = MatrixSpace(ZZ, m, n)
    _, basis = kernel(space(A))
    julia_form = Array(basis)
    return [ julia_form[:, j] for j in 1:size(julia_form, 2)]
end

"""
    lex_groebner_basis(A :: Matrix{Int})

Compute a lex Gröbner basis of `A` using the HNF algorithm.

It is assumed the toric ideal I_A is zero-dimensional.
TODO not sure this is the lex-gb I want. Review the theory on that!
"""
function lex_groebner_basis(
    A :: Matrix{Int}
)
    m, n = size(A)
    Mmn = MatrixSpace(ZZ, m, n)
    rnk, basis = kernel(Mmn(A))
    basis = basis' #Put the basis in row form
    hnf_basis = hnf(basis)
    normalize_hnf!(hnf_basis)
    return rnk, hnf_basis
end

"""
    truncate_markov(markov :: Vector{Vector{Int}}, instance :: IPInstance, truncation_type :: Symbol) :: Vector{Vector{Int}}

Compute a subset of `markov` including only vectors which shouldn't be truncated
according to the rule given by `truncation_type`.
"""
function truncate_markov(
    markov :: Vector{Vector{Int}},
    instance :: IPInstance,
    truncation_type :: Symbol
) :: Vector{Vector{Int}}
    truncated_basis = Vector{Int}[]
    should_truncate = truncation_type != :None
    for v in markov
        truncated = GBElements.truncate(
            v, instance.A, instance.b, instance.u,
            instance.model, instance.model_cons,
            should_truncate, truncation_type
        )
        if !truncated
            push!(truncated_basis, v)
        end
    end
    return truncated_basis
end

random_lifting(sigma) = rand(sigma)
minimum_lifting(sigma) = minimum(sigma)

"""
    project_and_lift(instance :: IPInstance, truncation_type :: Symbol) :: Vector{Vector{Int}}

Compute a Markov basis of `instance` using the project-and-lift algorithm.

Truncation is done with respect to `truncation_type`. If `truncation_type` is None, then
a full Markov basis is computed instead.

TODO I probably could make projections more efficient. But at least
this way they are very simple to implement...
"""
function project_and_lift(
    instance :: IPInstance;
    truncation_type :: Symbol = :None
) :: Vector{Vector{Int}}
    rank, hnf_basis = lex_groebner_basis(instance.A)
    @show hnf_basis
    #Now, the first `rank` columns of hnf_basis should be LI
    #and thus a Markov basis of the corresponding projection
    sigma = Vector(rank+1:instance.n)
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
    projection = nonnegativity_relaxation(instance, nonnegative)
    #Main loop: lift each variable in sigma
    while !isempty(sigma)
        i = minimum_lifting(sigma) #Pick some variable to lift
        perm_i = projection.permutation[i] #This is the index of i in projection
        if is_bounded(perm_i, projection)
            #A GB wrt the adequate order is a Markov basis of the lifted problem
            #Set up the right order in projection. minimize -i = maximize i
            update_objective!(projection, perm_i)
            alg = BuchbergerAlgorithm(
                markov, projection, truncation_type=truncation_type
            )
            markov = GBAlgorithms.run(alg, quiet=true)
        else
            #Find some vector u in ker(A) with u_i > 0 and u_{sigma_bar} >= 0
            u = unboundedness_proof(projection, nonnegative, perm_i)
            push!(markov, u)
        end
        #Finished lifting i, remove it from sigma
        i_index = findfirst(isequal(i), sigma)
        deleteat!(sigma, i_index)
        nonnegative[i] = true
        projection = nonnegativity_relaxation(instance, nonnegative)
        #Truncate the Markov basis
        markov = truncate_markov(markov, instance, truncation_type)
    end
    return markov
end

"""
    simple_markov(instance :: IPInstance) :: Vector{Vector{Int}}

Compute a Markov basis of `instance` with the simplified algorithm. Assumes all data in
instance.A and instance.b is non-negative.

The simplified algorithm uses the fact (proved in, e.g. Thomas and Weismantel (1997))
that the unit vectors on the original variables of the IP form a Markov basis.
"""
function simple_markov(
    instance::IPInstance
)::Vector{Vector{Int}}
    #Check whether the hypotheses hold
    @assert IPInstances.nonnegative_data_only(instance)
    #Build "trivial" Markov basis
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
Compute a Markov basis of an IP.
"""
function markov_basis end

function markov_basis(
    instance :: IPInstance;
    algorithm :: Symbol = :Any,
    truncation_type :: Symbol = :None
) :: Vector{Vector{Int}}
    if algorithm == :Any
        if IPInstances.nonnegative_data_only(instance)
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
        basis = project_and_lift(instance, truncation_type=truncation_type)
    end
    return basis
end

function markov_basis(
    A :: Array{Int, 2},
    b :: Vector{Int},
    C :: Array{Int, 2},
    u :: Vector{Int};
    algorithm :: Symbol = :Any,
    truncation_type :: Symbol = :None
) :: Vector{Vector{Int}}
    instance = IPInstance(A, b, C, u)
    return markov_basis(
        instance, algorithm=algorithm, truncation_type=truncation_type
    )
end

end
