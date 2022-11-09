"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis
"""
module Markov

export markov_basis, lex_groebner_basis

using AbstractAlgebra

using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.GBAlgorithms
using IPGBs.GBElements
using IPGBs.IPInstances
using IPGBs.Orders

"""
    normalize_hnf!(H :: Generic.MatSpaceElem{T})

Change `H` to an upper row HNF matrix satisfying the following property:
all entries above a pivot are non-positive and of smaller magnitude than the pivot

Assumes `H` is already in HNF as defined by AbstractAlgebra.jl, that is, it is in
upper row HNF satisfying:
- all entries above a pivot are non-negative and of smaller magnitude than the pivot
"""
function normalize_hnf!(
    H::Generic.MatSpaceElem{T}
) where {T}
    m, n = size(H)
    for i in 1:m
        #Find the pivot in row i
        j = i
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
    H::Generic.MatSpaceElem{T}
)::Bool where {T}
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
    group_relaxation_markov(instance::IPInstance)

Compute a Markov basis of the group relaxation of `instance`, returning this basis in
matrix form and its rank.

The computation is done using the Hermite Normal Form algorithm.
"""
function group_relaxation_markov(instance::IPInstance)
    basis = IPInstances.lattice_basis_projection(instance)
    uhnf_basis = hnf(basis)
    normalize_hnf!(uhnf_basis) #Normalize so that non-pivot entries are < 0
    return uhnf_basis, rank(uhnf_basis)
end

"""
    lex_groebner_basis(instance :: IPInstance) :: BinomialSet

Compute a lex Gröbner basis of `instance` using the HNF algorithm.

It is assumed the toric ideal of `instance.A` is zero-dimensional, i.e. that
the lattice optimization problem given by `instance` is of rank n.
"""
function lex_groebner_basis(
    instance::IPInstance
)::BinomialSet
    uhnf_basis, _ = group_relaxation_markov(instance)
    monomial_order = Orders.lex_order(instance.nonnegative_end)
    gb_vectors = Vector{Int}[]
    try
        #Convert lex GB elements to Vector{Int}
        gb_vectors = [
            reshape(Int.(Array(uhnf_basis[i, :])), ncols(uhnf_basis))
            for i in 1:nrows(uhnf_basis)
        ]
    catch e
        if isa(e, InexactError)
            error("Overflow in lex GB conversion.")
        else
            error("Unknown issue in converting lex GB to Vector{Int}.")
        end
    end
    lex_gb = BinomialSet(gb_vectors, monomial_order)
    return lex_gb
end

"""
    truncate_markov(markov :: Vector{Vector{Int}}, instance :: IPInstance, truncation_type :: Symbol) :: Vector{Vector{Int}}

Compute a subset of `markov` including only vectors which shouldn't be truncated
according to the rule given by `truncation_type`.
"""
function truncate_markov(
    markov::Vector{Vector{Int}},
    instance::IPInstance,
    truncation_type::Symbol
)::Vector{Vector{Int}}
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

minimum_lifting(sigma) = minimum(sigma)

struct ProjectAndLiftState
    instance::IPInstance #Instance being solved by P&L
    sigma::Vector{Int} #List of variables that still have to be lifted
    nonnegative::Vector{Bool} #Whether each variable is nonnegative at this state
    projection::IPInstance #Current projected lattice problem
    markov::Vector{Vector{Int}} #Current (partial) Markov basis
end

"""
    initialize_project_and_lift(instance :: IPInstance)

Create the initial state of the project-and-lift algorithm for IP `instance`, including
a Markov basis for its group relaxation.

The group relaxation Markov basis is obtained through the Hermite Normal Form algorithm.
"""
function initialize_project_and_lift(
    instance::IPInstance
)::ProjectAndLiftState
    hnf_basis, rnk = group_relaxation_markov(instance)
    #Now, the first `rank` columns of hnf_basis should be LI
    #and thus a Markov basis of the corresponding projection
    sigma = Vector(rnk+1:instance.n)
    nonnegative = nonnegative_vars(instance)
    for s in sigma
        nonnegative[s] = false
    end
    #This is a Markov basis of the projected problem
    markov = Vector{Int}[]
    for row in eachrow(Array(hnf_basis))
        v = Vector{Int}(row)
        push!(markov, lift_vector(v, instance))
    end
    projection = nonnegativity_relaxation(instance, nonnegative)
    @debug "Group relaxation Markov Basis: " markov
    @debug "Variables to lift: " sigma
    return ProjectAndLiftState(instance, sigma, nonnegative, projection, markov)
end

"""
    next(state :: ProjectAndLiftState; truncation_type :: Symbol) :: ProjectAndLiftState

Run a single iteration of the project-and-lift algorithm over `state`, returning the new
state after that iteration.

One iteration involves lifting a previously unlifted variable, either through a linear
program or a Gröbner Basis computation.
"""
function next(
    state::ProjectAndLiftState;
    truncation_type::Symbol = :None
)::ProjectAndLiftState
    i = minimum_lifting(state.sigma) #Pick some variable to lift
    perm_i = state.projection.permutation[i] #This is the index of i in projection
    u = unboundedness_proof(state.projection, state.nonnegative, perm_i)
    if isempty(u) #i is bounded in projection
        #Compute a GB in the adequate order
        @debug "Lifting $perm_i in bounded case, applying Buchberger's algorithm"
        update_objective!(state.projection, perm_i, state.sigma)
        alg = BuchbergerAlgorithm(
            state.markov, state.projection, truncation_type = truncation_type
        )
        markov = GBAlgorithms.run(alg, quiet = true)
        @debug "New Markov Basis obtained through Buchberger" markov
    else
        #u in ker(A) with u_i > 0 and u_{sigma_bar} >= 0
        @debug "Lifting $perm_i in unbounded case, add corresponding unbounded ray to the Markov Basis" u
        push!(markov, u)
    end
    #Finished lifting i, remove it from sigma
    i_index = findfirst(isequal(i), state.sigma)
    deleteat!(state.sigma, i_index)
    state.nonnegative[i] = true
    projection = nonnegativity_relaxation(state.instance, state.nonnegative)
    #Truncate the Markov basis
    markov = truncate_markov(markov, state.instance, truncation_type)
    return ProjectAndLiftState(
        state.instance, state.sigma, state.nonnegative, projection, markov
    )
end

"""
    is_finished(state::ProjectAndLiftState) :: Bool

Check whether the project-and-lift algorithm is ready to terminate at `state`.
"""
function is_finished(
    state::ProjectAndLiftState
)::Bool
    return isempty(state.sigma)
end

"""
    project_and_lift(instance :: IPInstance, truncation_type :: Symbol) :: Vector{Vector{Int}}

Compute a Markov basis of `instance` using the project-and-lift algorithm.

Truncation is done with respect to `truncation_type`. If `truncation_type` is None, then
a full Markov basis is computed instead.

TODO I probably could make projections more efficient. But at least
this way they are very simple to implement...
"""
function project_and_lift(
    instance::IPInstance;
    truncation_type::Symbol = :None
)::Vector{Vector{Int}}
    @debug "Initializing Project-and-lift"
    state = initialize_project_and_lift(instance)
    while !is_finished(state)
        l = length(state.sigma)
        @debug "Starting iteration with $l variables left to lift: " state.sigma
        state = next(state, truncation_type = truncation_type)
    end
    @debug "Ending P&L, found Markov Basis of length $(length(state.markov))"
    return state.markov
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
    @debug "Computing Markov basis using the Simple Markov algorithm"
    #Check whether the hypotheses hold
    @assert IPInstances.nonnegative_data_only(instance)
    #Build "trivial" Markov basis
    #Assumes there are upper bound constraints on the variables
    has_bounds = IPInstances.has_variable_bound_constraints(instance)
    n = size(instance.A, 2) - size(instance.A, 1) #Do not count slacks
    if has_bounds
        m = size(instance.A, 1) - n #Only count constraints that don't correspond to variable upper bounds
    else
        m = size(instance.A, 1)
    end
    basis = Vector{Int}[]
    g = Int[]
    for i in 1:n
        v = zeros(Int, n)
        v[i] = 1
        r = zeros(Int, n)
        r[i] = -1
        s = -copy(instance.A[1:m, i])
        if has_bounds
            g = vcat(v, s, r)
        else
            g = vcat(v, s)
        end
        push!(basis, g)
    end
    return basis
end

"""
Compute a Markov basis of an IP.
"""
function markov_basis end

function markov_basis(
    instance::IPInstance;
    algorithm::Symbol = :Any,
    truncation_type::Symbol = :None
)::Vector{Vector{Int}}
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
        basis = project_and_lift(instance, truncation_type = truncation_type)
    end
    return basis
end

function markov_basis(
    A::Array{Int,2},
    b::Vector{Int},
    C::Array{Int,2},
    u::Vector{Int};
    algorithm::Symbol = :Any,
    truncation_type::Symbol = :None
)::Vector{Vector{Int}}
    instance = IPInstance(A, b, C, u)
    return markov_basis(
        instance, algorithm = algorithm, truncation_type = truncation_type
    )
end

end
