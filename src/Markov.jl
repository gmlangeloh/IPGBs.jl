"""
Implementation of algorithms to compute Markov bases. These are necessary
as sets of generators of an ideal for computing Gröbner bases.

Two methods are implemented:
- a simplified method that only works when all data defining an IP is negative
- the project-and-lift algorithm of 4ti2 / Malkin's thesis
"""
module Markov

export markov_basis, project_and_lift

using IPGBs
using IPGBs.Binomials
using IPGBs.BinomialSets
using IPGBs.Buchberger
using IPGBs.GBAlgorithms
using IPGBs.GBElements
using IPGBs.GBTools
using IPGBs.IPInstances
using IPGBs.MatrixTools
using IPGBs.Orders
using IPGBs.SolverTools

using IPGBs.FourTi2

using JuMP

"""
    truncate_markov(
    markov::Vector{Vector{Int}},
    instance::IPInstance,
    truncation_type::Symbol
)::Vector{Vector{Int}}

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
        trunc_var_type, optimizer=instance.optimizer
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

"""
    most_negative(sigma, solution)

    Find the index of the most negative variable in `solution` indexed by `sigma`.
"""
function most_negative(sigma, solution)
    min_value = 0
    min_index = 0
    for i in sigma
        if solution[i] < min_value
            min_value = solution[i]
            min_index = i
        end
    end
    @assert min_value < 0
    return min_index
end

#
# ---------------------------------------
#

mutable struct ProjectAndLiftStats
    lift_time :: Float64
    optimize_time :: Float64
    num_lifts :: Int

    function ProjectAndLiftStats()
        new(0.0, 0.0, 0)
    end
end

function Base.show(io :: IO, stats :: ProjectAndLiftStats)
    println(io, "Project-and-Lift Stats: ")
    println(io, "lifts => ", stats.num_lifts)
    println(io, "lift_time => ", stats.lift_time)
    print(io, "optimize_time => ", stats.optimize_time)
end

function stats_after_lift(stats :: ProjectAndLiftStats, t :: Float64)
    stats.lift_time += t
    stats.num_lifts += 1
end

update_optimize_stats(stats :: ProjectAndLiftStats, t :: Float64) = stats.optimize_time += t

"""
    State of the project-and-lift algorithm representing the relevant
    information during a specific iteration.

    - `original_instance`: the original IPInstance being solved

    - `working_instance`: instance including potential upper bounds for efficiency

    - `unlifted`: the remaining unlifted variables (original indexing, following
    `working_instance`)

    - `nonnegative`: the variables that have nonnegative constraints (original indexing)

    - `relaxation`: IPInstance representing the problem with the non-negativity of
    the unlifted variables relaxed

    - `markov`: current partial Markov basis (indexed by the permuted variables)

    - `primal_solutions`: feasible solutions for the original problem (and thus all of
    its group relaxations). Indexed by the relaxation variables.

    - `dual_solution: feasible solution for `relaxation`. Computing a feasible
    solution is essentially free in project-and-lift, so we always compute one.
    If it is ever feasible for the original problem, it is also optimal for it.
    Uses the variable ordering of `relaxation`.

    - `optimal_solution`: The optimal solution for the original problem, if known
    (original indexing)

    - `has_optimal_solution`: Whether the optimal solution is known.

    - `stats`: Total statistics for this run of project-and-lift

    - `completion`: The GB algorithm to be used. :IPGBs and :4ti2 are supported.
"""
mutable struct ProjectAndLiftState
    original_instance::IPInstance
    working_instance::IPInstance
    unlifted::Vector{Int}
    nonnegative::Vector{Bool}
    relaxation::IPInstance
    markov::Vector{Vector{Int}}
    primal_solutions::Vector{Vector{Int}}
    dual_solution::Vector{Int}
    optimal_solution::Vector{Int}
    has_optimal_solution::Bool
    stats::ProjectAndLiftStats
    completion::Symbol
    optimize::Bool
    truncation_type::Symbol
    quiet::Bool
end

function ProjectAndLiftState(
    previous_state :: ProjectAndLiftState,
    new_relaxation :: IPInstance,
    new_markov :: Vector{Vector{Int}},
    new_primals :: Vector{Vector{Int}},
    new_dual :: Vector{Int}
)
    return ProjectAndLiftState(
        previous_state.original_instance, previous_state.working_instance,
        previous_state.unlifted, previous_state.nonnegative, new_relaxation,
        new_markov, new_primals, new_dual, previous_state.optimal_solution,
        previous_state.has_optimal_solution, previous_state.stats,
        previous_state.completion, previous_state.optimize,
        previous_state.truncation_type, previous_state.quiet
    )
end

function set_optimal_solution!(
    pl_state :: ProjectAndLiftState,
    solution :: Vector{Int}
)
    pl_state.optimal_solution = solution
    pl_state.has_optimal_solution = true
end

function Base.show(io :: IO, state :: ProjectAndLiftState)
    print(io, "Project-and-Lift State for instance: ", state.original_instance, "\n")
    print(io, "Unlifted = ", state.unlifted, "\n")
    print(io, "Markov = ", state.markov, "\n")
end

print_stats(s :: ProjectAndLiftState, quiet = false) = if !quiet println(s.stats) end

#
# Lifting feasible solutions using a Markov basis
#

function initial_solution(
    instance :: IPInstance,
    markov :: Vector{Vector{Int}}
) :: Vector{Int}
    solution = copy(fiber_solution(instance))
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
    #TODO: Refactor, use optimize_with!
    old_obj = instance.C[1, :]
    obj = SolverTools.cone_element(markov, optimizer=instance.optimizer)
    instance.C[1, :] = obj
    bs = BinomialSet(markov, instance, Binomial)
    binomial_sol = to_gbelement(solution, bs.order, Binomial, false)
    initialize_binomials(instance, bs.order)
    BinomialSets.reduce!(binomial_sol, bs, is_monomial_reduction=true)
    copyto!(solution, element(binomial_sol))
    instance.C[1, :] = old_obj
    @assert is_feasible_solution(instance, solution)
    return solution
end

"""
    relaxation_feasible(pl_state :: ProjectAndLiftState)

    The feasible solution for the relaxation, in the original variable ordering.
"""
function relaxation_feasible(pl_state :: ProjectAndLiftState)
    return pl_state.dual_solution[pl_state.relaxation.inverse_permutation]
end

function lift_solution_unbounded!(
    solution :: Vector{Int},
    ray :: Vector{Int}
)
    k = maximum([ceil(Int, -solution[i] / ray[i])
        for i in eachindex(ray) if ray[i] > 0 && solution[i] < 0],
        init = 0
    )
    solution .= solution .+ k * ray
end

#
# ----------------------------------------------
#

project_to_instance(v :: Vector{Int}, instance :: IPInstance) = v[1:instance.n]

"""
    optimize_relaxation!(pl_state :: ProjectAndLiftState)

    Optimize `pl_state.dual_solution`, which is feasible for `pl_state.relaxation`,
    obtaining an optimal solution for the relaxation.

    Gröbner bases are used to optimize the solution.
"""
function optimize_relaxation!(pl_state :: ProjectAndLiftState)
    if !pl_state.optimize || pl_state.has_optimal_solution
        return
    end
    is_optimal = false
    opt_solution = pl_state.optimal_solution
    all_solutions = Vector{Int}[]
    #TODO: Question: Should I actually use the primal solutions here?
    # They are obviously feasible for the relaxation, but I'm not sure
    # optimizing them here does anything...
    # - in particular, they may not be feasible for the original problem
    for solution in pl_state.primal_solutions
        push!(all_solutions, solution)
    end
    push!(all_solutions, pl_state.dual_solution)
    #TODO: For now, I won't try to reuse this Gröbner basis. Might be worthwhile
    # to try doing that eventually, though (easier Markov bases later?)
    #TODO: The Binary Truncation Criterion as currently implemented cannot be used
    # in group relaxations. Maybe it can only be used for variables with
    # non-negative constraints? This is an interesting extension to explore.
    # quick_truncation appears to have a similar problem...
    _, t, _, _, _ = @timed groebner_basis(
        pl_state.relaxation, pl_state.markov, solutions=all_solutions,
        quiet=pl_state.quiet, truncation_type=pl_state.truncation_type,
        use_binary_truncation=false, use_quick_truncation=false
    )
    update_optimize_stats(pl_state.stats, t)
    orig_dual = pl_state.dual_solution[pl_state.relaxation.inverse_permutation]
    # Any dual solution feasible for the original problem is optimal.
    if is_feasible_solution(pl_state.working_instance, orig_dual)
        opt_solution = project_to_instance(orig_dual, pl_state.original_instance)
        is_optimal = true
    end
    # TODO: Keep truncation bounds updated in opt_instance and relaxation
    # Need to update both opt_instance and corresponding solution slacks...
    #bkv = minimum([opt_instance.C[1, :]' * s for s in primal_solutions])
    if is_optimal
        set_optimal_solution!(pl_state, opt_solution)
    end
end

"""
    initialize_project_and_lift(instance :: IPInstance)

Create the initial state of the project-and-lift algorithm for IP `instance`, including
a Markov basis for its group relaxation.

The group relaxation Markov basis is obtained through the Hermite Normal Form algorithm.
"""
function initialize_project_and_lift(
    instance :: IPInstance;
    completion :: Symbol = :IPGBs,
    optimize :: Bool = false,
    truncation_type :: Symbol = :LP,
    quiet :: Bool = true,
    solution :: Vector{Int} = zeros(Int, instance.n)
)::ProjectAndLiftState
    # Update instance to include upper bounds for efficiency if possible
    opt_instance = instance
    primal_solutions = Vector{Int}[]
    #TODO: Add upper bounds later!!! Currently, this doesn't work
    #if optimize && is_feasible_solution(instance, solution)
    #    # In this case, we should use the given solution as an upper bound
    #    # GBs are smaller this way.
    #    val = round(Int, instance.C[1, :]' * solution)
    #    opt_instance = add_constraint(instance, round.(Int, instance.C[1, :]), val)
    #    push!(primal_solutions, solution)
    #end
    # Compute a group relaxation with its corresponding Markov basis
    init_basis_algorithm = optimize ? :SimplexBasis : :Any
    uhnf_basis, proj_basis, unlifted = lattice_basis_projection(
        opt_instance, init_basis_algorithm
    )
    nonnegative = nonnegative_vars(opt_instance)
    for s in unlifted
        nonnegative[s] = false
    end
    markov = Vector{Int}[]
    for row in eachrow(Array(uhnf_basis))
        v = Vector{Int}(row)
        lifted = lift_vector(v, proj_basis, lattice_basis(opt_instance))
        push!(markov, lifted)
    end
    relaxation = nonnegativity_relaxation(opt_instance, nonnegative)
    permuted_markov = apply_permutation(markov, relaxation.permutation)
    primal_solutions = apply_permutation(primal_solutions, relaxation.permutation)
    #Find initial dual solution if possible
    relaxation_solution = initial_solution(relaxation, permuted_markov)
    pl_state = ProjectAndLiftState(
        instance, opt_instance, unlifted, nonnegative, relaxation,
        permuted_markov, primal_solutions, relaxation_solution,
        zeros(Int, instance.n), false, ProjectAndLiftStats(),
        completion, optimize, truncation_type, quiet
    )
    optimize_relaxation!(pl_state)
    return pl_state
end

"""
    relaxation_index(i :: Int, relaxation :: IPInstance)

    Index of original variable `i` in `relaxation`.
"""
relaxation_index(i :: Int, r :: IPInstance) = r.inverse_permutation[i]
relaxation_index(v :: Vector{Int}, r :: IPInstance) = [ relaxation_index(i, r) for i in v]

"""
    original_index(i :: Int, r :: IPInstance)

    Index of `r`'s variable `i` in the original instance.
"""
original_index(i :: Int, r :: IPInstance) = r.permutation[i]

"""
    can_lift(markov :: Vector{Vector{Int}}, i :: Int)

    Return true if the variable indexed by `i` can be immediately lifted given
    the current Markov basis `markov`.

    This happens when there is no element in `markov` with a positive coefficient.
"""
function can_lift(markov :: Vector{Vector{Int}}, i :: Int)
    return all(v[i] <= 0 for v in markov)
end

"""
    lift_variables!(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
)

    Mark all unlifted variables that can be lifted in `pl_state` as nonnegative
    and lift them.
"""
function lift_variables!(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
)
    i = 1
    while i <= length(pl_state.unlifted)
        variable = pl_state.unlifted[i]
        #Markov is permuted here, so we need to take that into account
        rel_idx = relaxation_index(variable, pl_state.relaxation)
        if can_lift(markov, rel_idx)
            pl_state.nonnegative[variable] = true
            deleteat!(pl_state.unlifted, i)
        else
            i += 1
        end
    end
end

"""
    lift_bounded(
    pl_state :: ProjectAndLiftState,
    i :: Int
) :: Vector{Vector{Int}}

    Lift the variable indexed by `i` in `pl_state` using a Gröbner basis computation.
    This generates a Markov basis for the next group relaxation.
"""
function lift_bounded(
    pl_state :: ProjectAndLiftState,
    i :: Int
) :: Vector{Vector{Int}}
    #Compute a GB in the adequate monomial order
    @info "P&L bounded case" i length(pl_state.markov)
    update_objective!(
        pl_state.relaxation, i, relaxation_index(pl_state.unlifted, pl_state.relaxation)
    )
    if pl_state.completion == :IPGBs
        #This will automatically lift pl.dual_solutions
        #TODO: Using quick truncation breaks something here. It may be a bug in my code,
        #or some implicit hypothesis that I'm not meeting. Figure this out later!
        markov = IPGBs.groebner_basis(
            pl_state.relaxation, pl_state.markov, solutions=[pl_state.dual_solution],
            truncation_type=pl_state.truncation_type, use_quick_truncation=false,
            quiet = pl_state.quiet, use_binary_truncation=false
        )
    elseif pl_state.completion == :FourTi2
        gb, sol, _ = FourTi2.groebnernf(
            pl_state.relaxation, pl_state.markov, pl_state.dual_solution
        )
        markov = GBTools.tovector(gb)
        copyto!(pl_state.dual_solution, sol)
    else
        error("Unknown completion algorithm: $(pl_state.completion)")
    end
    @info "Markov size after lifting: " length(markov)
    return markov
end

"""
    lift_unbounded(
    pl_state :: ProjectAndLiftState,
    i :: Int,
    ray :: Vector{Int}
) :: Vector{Vector{Int}}

    Lift the variable indexed by `i` in `pl_state` using `ray` that proves that
    the variable is unbounded in `pl_state.relaxation`.
"""
function lift_unbounded(
    pl_state :: ProjectAndLiftState,
    i :: Int,
    ray :: Vector{Int}
) :: Vector{Vector{Int}}
    @info "P&L unbounded case" ray
    if !(ray in pl_state.markov)
        push!(pl_state.markov, ray)
        lift_solution_unbounded!(pl_state.dual_solution, ray)
    end
    i_index = findfirst(isequal(i), pl_state.unlifted)
    deleteat!(pl_state.unlifted, i_index)
    pl_state.nonnegative[i] = true
    return pl_state.markov
end

"""
    relax_and_reorder(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
) :: Tuple{IPInstance, Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Int}}

    Relax the non-negativity constraints of `pl_state.working_instance` according to
    `pl_state.nonnegative`, and reorder the Markov basis and solutions to match the
    new variable ordering.
"""
function relax_and_reorder(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
) :: Tuple{IPInstance, Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Int}}
    #Before taking the relaxation, we bring the Markov basis and known solutions back into
    #the original variable ordering
    lifted_markov = apply_permutation(markov, pl_state.relaxation.inverse_permutation)
    dual_solution = pl_state.dual_solution[pl_state.relaxation.inverse_permutation]
    primal_solutions = apply_permutation(
        pl_state.primal_solutions, pl_state.relaxation.inverse_permutation
    )

    #Now we take the relaxation wrt the non-negativity constraints specified in pl_state
    relaxation = nonnegativity_relaxation(pl_state.working_instance, pl_state.nonnegative)

    #And we permute the Markov basis and solutions back to the new variable ordering
    lifted_markov = apply_permutation(lifted_markov, relaxation.permutation)
    dual_solution = dual_solution[relaxation.permutation]
    primal_solutions = apply_permutation(pl_state.primal_solutions, relaxation.permutation)
    return relaxation, lifted_markov, primal_solutions, dual_solution
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
    pl_state :: ProjectAndLiftState
)::ProjectAndLiftState
    #Find out which variable to lift next
    if pl_state.optimize
        #Lift the variable with the most negative value in the dual solution
        #Heuristically, this should help making the dual solution feasible,
        #leading to early termination
        relaxation_unlifted = relaxation_index(pl_state.unlifted, pl_state.relaxation)
        relaxation_i = most_negative(relaxation_unlifted, pl_state.dual_solution)
        i = original_index(relaxation_i, pl_state.relaxation)
    else
        #Any variable is fine, lift the first one
        i = first_variable(pl_state.unlifted)
        relaxation_i = relaxation_index(i, pl_state.relaxation)
    end
    #Now figure out whether i is bounded or unbounded by linear programming
    ray = unboundedness_proof(pl_state.relaxation, relaxation_i)
    if isempty(ray) #No ray means the variable is bounded
        lifted_markov, t, _, _, _ = @timed lift_bounded(pl_state, relaxation_i)
    else
        lifted_markov, t, _, _, _ = @timed lift_unbounded(pl_state, i, ray)
    end
    stats_after_lift(pl_state.stats, t)
    return lift_and_relax(pl_state, lifted_markov)
end

"""
    lift_and_relax(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
) :: ProjectAndLiftState

    Lift the variables in `pl_state` according to `markov` generating a new state
    with the updated group relaxation, Markov basis and solutions.
"""
function lift_and_relax(
    pl_state :: ProjectAndLiftState,
    markov :: Vector{Vector{Int}}
) :: ProjectAndLiftState
    lift_variables!(pl_state, markov)
    relaxation, lifted_markov, primals, dual = relax_and_reorder(pl_state, markov)
    lifted_markov = truncate_markov(lifted_markov, relaxation, pl_state.truncation_type)
    new_state = ProjectAndLiftState(pl_state, relaxation, lifted_markov, primals, dual)
    optimize_relaxation!(new_state)
    return new_state
end

all_lifted(pl_state :: ProjectAndLiftState) = isempty(pl_state.unlifted)

"""
    is_finished(pl_state::ProjectAndLiftState)

Check whether the project-and-lift algorithm is ready to terminate at `pl_state`.
"""
function is_finished(
    pl_state::ProjectAndLiftState
)
    return all_lifted(pl_state) || pl_state.has_optimal_solution
end

"""
    project_and_lift(
    instance::IPInstance;
    completion :: Symbol = :IPGBs,
    truncation_type::Symbol = :LP,
    optimize :: Bool = false,
    feasible :: Bool = false,
    quiet :: Bool = true,
    solution :: Vector{Int} = zeros(Int, instance.n)
)::Tuple{Vector{Vector{Int}}, Bool, Vector{Int}, Float64}

Compute a Markov basis of `instance` using the project-and-lift algorithm.

Truncation is done with respect to `truncation_type`. If `truncation_type` is None, then
a full Markov basis is computed instead.
"""
function project_and_lift(
    instance::IPInstance;
    completion :: Symbol = :IPGBs,
    truncation_type::Symbol = :LP,
    optimize :: Bool = false,
    feasible :: Bool = false,
    quiet :: Bool = true,
    solution :: Vector{Int} = zeros(Int, instance.n)
)::Tuple{Vector{Vector{Int}}, Bool, Vector{Int}, Float64}
    #Compute the initial group relaxation and lift as many variables as possible
    pl_state = initialize_project_and_lift(
        instance, completion=completion, truncation_type=truncation_type,
        optimize=optimize, quiet=quiet, solution=solution
    )
    pl_state = lift_and_relax(pl_state, pl_state.markov)

    #Main loop: lift all remaining variables via LP or GBs
    # Since computing a feasible solution is essentially free, we do that as well
    while !is_finished(pl_state)
        pl_state = next(pl_state)
    end

    #Postprocessing: check solutions are consistent, compute optimal value
    @assert all(is_feasible_solution(instance, solution, pl_state.relaxation.inverse_permutation) for solution in pl_state.primal_solutions)
    #@assert is_feasible_solution(instance, relaxation_feasible(pl_state))
    #@assert !pl_state.has_optimal_solution || is_feasible_solution(instance, pl_state.optimal_solution)
    if feasible && !pl_state.has_optimal_solution
        pl_state.optimal_solution = relaxation_feasible(pl_state)
        pl_state.has_optimal_solution = true
    end
    print_stats(pl_state, quiet)
    optimal_value = instance.C[1, :]' * pl_state.optimal_solution
    return pl_state.markov, pl_state.has_optimal_solution, pl_state.optimal_solution, optimal_value
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
        #Only count constraints that don't correspond to variable upper bounds
        m = size(instance.A, 1) - n
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
    solution :: Vector{Int} = zeros(Int, instance.n),
    algorithm::Symbol = :Any,
    truncation_type::Symbol = :LP,
    optimize::Bool = false,
    quiet::Bool = true
)::Vector{Vector{Int}}
    @debug "Starting to compute Markov basis for " instance
    if !is_bounded(instance)
        return [zeros(Int, instance.n)]
    end
    can_use_simple = nonnegative_data_only(instance) && has_slacks(instance)
    has_solution = !iszero(solution) && is_feasible_solution(instance, solution)
    if algorithm == :Any
        if can_use_simple && !has_solution
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
        basis, _, _, _ = project_and_lift(
            instance, solution=solution, truncation_type=truncation_type,
            optimize=optimize, quiet=quiet
        )
    end
    return basis
end

function markov_basis(
    A::Array{Int,2},
    b::Vector{Int},
    C::Array{Int,2},
    u::Vector{Int};
    algorithm::Symbol = :Any,
    truncation_type::Symbol = :None,
    quiet :: Bool = true
)::Vector{Vector{Int}}
    instance = IPInstance(A, b, C, u)
    return markov_basis(
        instance, algorithm = algorithm, truncation_type = truncation_type, quiet = quiet
    )
end

function markov_basis(
    model :: JuMP.Model;
    algorithm :: Symbol = :Any,
    truncation_type :: Symbol = :None,
    quiet :: Bool = true
)
    instance = IPInstance(model)
    return markov_basis(
        instance, algorithm = algorithm, truncation_type = truncation_type, quiet = quiet
    )
end

function markov_basis(
    filename :: String;
    algorithm :: Symbol = :Any,
    truncation_type :: Symbol = :None,
    quiet :: Bool = true
)
    instance = IPInstance(filename)
    return markov_basis(
        instance, algorithm = algorithm, truncation_type = truncation_type, quiet = quiet
    )
end

end
