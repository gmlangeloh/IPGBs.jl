"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis
"""
module Markov

export markov_basis

using MIPMatrixTools.GBTools
using MIPMatrixTools.IPInstances
using MIPMatrixTools.MatrixTools
using MIPMatrixTools.SolverTools

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.GBAlgorithms
using IPGBs.GBElements
using IPGBs.Orders

using IPGBs.FourTi2

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
    if !should_truncate
        return markov
    end
    trunc_var_type = Int
    if truncation_type == :LP
        trunc_var_type = Real
    end
    model, _, constrs = SolverTools.feasibility_model(
        instance.A, instance.b, instance.u, nonnegative_vars(instance),
        trunc_var_type
    )
    for v in markov
        truncated = GBElements.truncate(
            v, instance.A, instance.b, instance.u,
            model, constrs, should_truncate, truncation_type
        )
        if !truncated
            push!(truncated_basis, v)
        end
    end
    return truncated_basis
end

#
# How to pick the next variable to lift?
#

first_variable(sigma) = minimum(sigma)
most_negative(sigma, solution) = argmin(solution[sigma])

#
# ---------------------------------------
#

"""
    State of the project-and-lift algorithm representing the relevant
    information during a specific iteration.
    `original_instance`: the original IPInstance being solved
    `working_instance`: instance including potential upper bounds for efficiency
    `unlifted`: the remaining unlifted variables (original indexing, following `working_instance`)
    `nonnegative`: the variables that have nonnegative constraints (original indexing)
    `relaxation`: IPInstance representing the problem with the non-negativity of
    the unlifted variables relaxed
    `markov`: current partial Markov basis (indexed by the permuted variables)
    `primal_solutions`: feasible solutions for the original problem (and thus all of
    its group relaxations)
    `dual_solutions`: feasible solutions for `relaxation`. Computing a feasible
    solution is essentially free in project-and-lift, so we always compute one.
    If it is ever feasible for the original problem, it is also optimal for it.
    `optimal_solution`: The optimal solution for the original problem, if known.
"""
struct ProjectAndLiftState
    original_instance::IPInstance
    working_instance::IPInstance
    unlifted::Vector{Int}
    nonnegative::Vector{Bool}
    relaxation::IPInstance
    markov::Vector{Vector{Int}}
    primal_solutions::Vector{Vector{Int}}
    dual_solutions::Vector{Vector{Int}}
    optimal_solution::Vector{Int}
end

function Base.show(io :: IO, state :: ProjectAndLiftState)
    print(io, "Project-and-Lift State for instance: ", state.original_instance, "\n")
    print(io, "Unlifted = ", state.unlifted, "\n")
    print(io, "Markov = ", state.markov, "\n")
end

#
# Lifting feasible solutions using a Markov basis
#

function initial_solution(
    instance :: IPInstance,
    markov :: Vector{Vector{Int}}
) :: Vector{Int}
    solution = copy(instance.fiber_solution)
    #Make the solution feasible for the initial state (group relaxation)
    for i in eachindex(solution)
        if IPInstances.is_nonnegative(i, instance) && solution[i] < 0
            #Find Markov basis element with positive coefficient on i
            #And add them to solution as many times as needed to guarantee
            #solution[i] >= 0
            for m in markov
                if m[i] > 0
                    k = ceil(Int, -solution[i] / m[i])
                    solution += k * m
                    break
                end
            end
        end
    end
    @assert is_feasible_solution(instance, solution)
    return solution
end

function lift_solution_unbounded(
    solution :: Vector{Int}, 
    ray :: Vector{Int}
) :: Vector{Int}
    k = maximum([ceil(Int, -solution[i] / ray[i]) 
        for i in eachindex(ray) if ray[i] > 0 && solution[i] < 0]
    )
    return solution + k * ray
end

#
# ----------------------------------------------
#

project_to_instance(v :: Vector{Int}, instance :: IPInstance) = v[1:instance.n]

function optimize_with_markov(
    opt_instance :: IPInstance,
    primal_solutions :: Vector{Vector{Int}},
    dual_solutions :: Vector{Vector{Int}},
    relaxation :: IPInstance,
    permuted_markov :: Vector{Vector{Int}};
    optimize :: Bool = false
)
    is_optimal = false
    opt_solution = zeros(Int, opt_instance.n)
    if optimize
        all_solutions = [ primal_solutions; dual_solutions ]
        # For now, I won't try to reuse this Gröbner basis. Might be worthwhile
        # to try doing that eventually, though (easier Markov bases later?)
        _ = groebner_basis(permuted_markov, relaxation, all_solutions)
        for relaxation_s in dual_solutions
            s = relaxation_s[relaxation.inverse_permutation]
            # Any dual solution feasible for the original problem is optimal.
            if is_feasible_solution(opt_instance, s)
                opt_solution = s
                is_optimal = true
                break
            end
        end
        #TODO: Keep truncation bounds updated in opt_instance and relaxation
        bkv = minimum([opt_instance.C[1, :]' * s for s in primal_solutions])
        #Need to update both opt_instance and corresponding solution slacks...
    end
    return opt_solution, is_optimal
end

"""
    initialize_project_and_lift(instance :: IPInstance)

Create the initial state of the project-and-lift algorithm for IP `instance`, including
a Markov basis for its group relaxation.

The group relaxation Markov basis is obtained through the Hermite Normal Form algorithm.
"""
function initialize_project_and_lift(
    instance :: IPInstance;
    optimize :: Bool = false,
    solution :: Vector{Int} = zeros(Int, instance.n),
    best_value :: Ref{Int} = Ref(0)
)::ProjectAndLiftState
    # Update instance to include upper bounds for efficiency if possible
    opt_instance = instance
    if optimize && is_feasible_solution(instance, solution)
        # In this case, we should use the given solution as an upper bound
        # GBs are smaller this way.
        opt_instance = add_constraint(instance, instance.C[1, :], best_value[] - 1)
    end
    # Compute a group relaxation with its corresponding Markov basis
    basis, uhnf_basis, unlifted = lattice_basis_projection(opt_instance)
    nonnegative = nonnegative_vars(opt_instance)
    for s in unlifted
        nonnegative[s] = false
    end
    markov = [ lift_vector(v, basis, opt_instance) 
        for v in eachrow(Array(uhnf_basis))
    ]
    relaxation = nonnegativity_relaxation(opt_instance, nonnegative)
    permuted_markov = IPInstances.apply_permutation(markov, relaxation.permutation)
    #Find initial primal and dual solutions if possible
    relaxation_solution = initial_solution(relaxation, permuted_markov)
    primal_solutions = Vector{Int}[]
    if !iszero(solution)
        push!(primal_solutions, solution)
    end
    dual_solutions = [ relaxation_solution ]
    optimize_with_markov()
    return ProjectAndLiftState(
        instance, opt_instance, unlifted, nonnegative, relaxation, 
        permuted_markov, primal_solutions, dual_solutions
    )
end

"""
    relaxation_index(i :: Int, relaxation :: IPInstance)

    Index of original variable `i` in `relaxation`.
"""
function relaxation_index(i :: Int, relaxation :: IPInstance)
    return relaxation.inverse_permutation[i]
end

"""
    can_lift(markov :: Vector{Vector{Int}}, i :: Int)

    Return true if the variable indexed by `i` can be immediately lifted given
    the current Markov basis `markov`.

    This happens when there is no element in `markov` with a positive coefficient.
"""
function can_lift(markov :: Vector{Vector{Int}}, i :: Int)
    return all(v[i] <= 0 for v in markov)
end

function lift_variables!(
    pl :: ProjectAndLiftState;
    markov :: Union{Vector{Vector{Int}}, Nothing} = nothing
)
    if isnothing(markov)
        markov = pl.markov
    end
    i = 1
    while i <= length(pl.unlifted)
        variable = pl.unlifted[i]
        #Markov is permuted here, so we need to take that into account
        if can_lift(markov, relaxation_index(variable, pl.relaxation))
            pl.nonnegative[variable] = true
            deleteat!(pl.unlifted, i)
        else
            i += 1
        end
    end
end

"""
    bounded_case(
    s :: ProjectAndLiftState, 
    i :: Int;
    completion :: Symbol = :Buchberger,
    truncation_type :: Symbol = :None
) :: Vector{Vector{Int}}


"""
function lift_bounded(
    pl :: ProjectAndLiftState, 
    i :: Int;
    completion :: Symbol = :Buchberger,
    truncation_type :: Symbol = :None
) :: Vector{Vector{Int}}
    #Compute a GB in the adequate monomial order
    @info "P&L bounded case" i length(pl.markov)
    update_objective!(pl.relaxation, i, relaxation_index.(pl.unlifted, pl.relaxation))
    if completion == :Buchberger
        #This will automatically lift s.solution
        markov = IPGBs.groebner_basis(pl.markov, pl.relaxation, [pl.solution], truncation_type=truncation_type)
    elseif completion == :FourTi2
        gb, sol, _ = FourTi2.groebnernf(pl.relaxation, pl.markov, pl.solution)
        markov = GBTools.tovector(gb)
        copyto!(pl.solution, sol)
    else
        error("Unknown completion algorithm: $completion")
    end
    lift_variables!(pl, markov=markov)
    @info "Markov size after lifting: " length(markov)
    return markov
end

"""
    unbounded_case(s :: ProjectAndLiftState, i :: Int, ray :: Vector{Int}) :: Vector{Vector{Int}}


"""
function lift_unbounded(s :: ProjectAndLiftState, i :: Int, ray :: Vector{Int}) :: Vector{Vector{Int}}
    @info "P&L unbounded case" ray
    if !(ray in s.markov)
        push!(s.markov, ray)
    end
    i_index = findfirst(isequal(i), s.unlifted)
    deleteat!(s.unlifted, i_index)
    s.nonnegative[i] = true
    return s.markov
end

"""
    next(
    s::ProjectAndLiftState;
    completion :: Symbol = :Buchberger,
    truncation_type::Symbol = :LP
)::ProjectAndLiftState

Run a single iteration of the project-and-lift algorithm over `state`, 
returning the new state after that iteration.

One iteration involves lifting a previously unlifted variable, either through 
a linear program or a Gröbner Basis computation.
"""
function next(
    pl::ProjectAndLiftState;
    completion :: Symbol = :Buchberger,
    truncation_type::Symbol = :LP
)::ProjectAndLiftState
    i = first_variable(pl.unlifted) #Pick some variable to lift
    relaxation_i = relaxation_index(i, pl.relaxation)
    ray = unboundedness_proof(pl.relaxation, relaxation_i)
    if isempty(ray)
        lifted_markov = lift_bounded(
            pl, relaxation_i, completion=completion,
            truncation_type=truncation_type
        )
        lifted_solution = pl.solution
    else
        lifted_markov = lift_unbounded(pl, i, ray)
        lifted_solution = lift_solution_unbounded(pl.solution, ray)
    end
    #TODO: Clean all this stuff up, it looks terrible!!!
    #markov needs to be permuted to match the order of variables in the new
    #group relaxation. To do this, we first put it back to the original variable
    #order and then apply the permutation of the relaxation
    lifted_markov = IPInstances.apply_permutation(lifted_markov, pl.relaxation.inverse_permutation)
    lifted_solution = lifted_solution[pl.relaxation.inverse_permutation]
    next_relaxation = nonnegativity_relaxation(pl.working_instance, pl.nonnegative)
    lifted_markov = IPInstances.apply_permutation(lifted_markov, next_relaxation.permutation)
    lifted_solution = lifted_solution[next_relaxation.permutation]
    lifted_markov = truncate_markov(lifted_markov, next_relaxation, truncation_type)
    return ProjectAndLiftState(
        pl.original_instance, pl.working_instance, pl.unlifted, pl.nonnegative, 
        next_relaxation, lifted_markov, lifted_solution
    )
end

"""
    is_finished(state::ProjectAndLiftState) :: Bool

Check whether the project-and-lift algorithm is ready to terminate at `state`.
"""
function is_finished(
    state::ProjectAndLiftState
)::Bool
    return isempty(state.unlifted)
end

"""
    project_and_lift(instance :: IPInstance, truncation_type :: Symbol) :: Vector{Vector{Int}}

Compute a Markov basis of `instance` using the project-and-lift algorithm.

Truncation is done with respect to `truncation_type`. If `truncation_type` is None, then
a full Markov basis is computed instead.
"""
function project_and_lift(
    instance::IPInstance;
    completion :: Symbol = :Buchberger,
    truncation_type::Symbol = :LP,
    optimize :: Bool = false,
    solution :: Vector{Int} = zeros(Int, instance.n),
    best_value :: Ref{Int} = Ref(0)
)::Vector{Vector{Int}}
    @debug "Initializing Project-and-lift"
    pl = initialize_project_and_lift(
        instance, optimize=optimize, solution=solution, best_value=best_value
    )
    #Lift as many variables as possible before starting
    lift_variables!(pl)
    #Main loop: lift all remaining variables via LP or GBs
    while !is_finished(pl)
        @debug "Starting iteration with $(length(pl.sigma)) variables left to lift: "
        pl = next(pl, completion=completion, truncation_type=truncation_type)
    end
    @assert is_feasible_solution(instance, pl.solution)
    @debug "Ending P&L, found Markov Basis of length $(length(pl.markov))"
    #Update solution and best known value
    copyto!(solution, pl.solution)
    best_value[] = instance.C[1, :]' * solution
    return pl.markov
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
    @assert IPInstances.has_slacks(instance)
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
    truncation_type::Symbol = :LP
)::Vector{Vector{Int}}
    @debug "Starting to compute Markov basis for " instance
    if algorithm == :Any
        if IPInstances.nonnegative_data_only(instance) && IPInstances.has_slacks(instance)
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

function optimize(
    instance :: IPInstance;
    completion :: Symbol = :Buchberger,
    truncation_type :: Symbol = :LP
) :: Tuple{Vector{Int}, Int}
    solution = copy(instance.fiber_solution)
    value = instance.C[1, :]' * solution
    project_and_lift(
        instance, completion=completion, truncation_type=truncation_type,
        optimize=true, solution=solution, best_value=Ref(value)
    )
    return solution, value
end

end
