"""
A module to solve multi-objective integer programming problems using
epsilon-constraint methods and test sets / Gröbner Bases.
"""
module MultiObjectiveAlgorithms
export moip_gb_solve, moip_walkback

using IPGBs
using IPGBs.BinomialSets
using IPGBs.FourTi2
using IPGBs.IPInstances
using IPGBs.MatrixTools
using IPGBs.Orders
using IPGBs.Walkback

using IPGBs.MultiObjectiveTools
using IPGBs.SingleObjective
using IPGBs.MultiObjectiveStats

using JuMP
using GLPK
using LinearAlgebra
using Random

"""
    new_step_list(steps :: Vector{Vector{Int}}, slacks :: Vector{Int}) :: Vector{Vector{Int}}

Add the `slacks` vector to the `steps` list if it is not dominated by any vector
in `steps`, and remove from `steps` any vector dominated by `slacks`.

In the biobjective case this is unnecessary, but dropping useless steps
can lead to better performance when there are 3+ objectives.
"""
function new_step_list(
    steps::Vector{Vector{Int}},
    slacks::Vector{Int}
)::Vector{Vector{Int}}
    old_steps = steps
    #Remove everything in steps that is dominated by slacks
    steps = filter(step -> !all(step .>= slacks), steps)
    for step in steps
        if all(step .<= slacks) #slacks is dominated by step, ...
            #...it doesn't have to be added to steps
            return old_steps
        end
    end
    push!(steps, copy(slacks))
    return steps
end

"""
    minimalsteps(
    gb::Vector{Vector{Int}},
    x::Vector{Int}
)::Vector{Vector{Int}}

Return the set of minimal steps in objective value from `x` that may yield new
points of the Pareto frontier in an epsilon-constraint method.
"""
function minimal_steps(
    gb::Vector{Vector{Int}},
    x::Vector{Int},
    variable_indices::Vector{Int},
    slack_indices::Vector{Int}
)::Vector{Vector{Int}}
    if isempty(gb)
        return Vector{Int}[]
    end
    num_slacks = length(slack_indices)
    first_slack_idx = slack_indices[1]
    steps = [fill(typemax(Int64), num_slacks)]
    slacks = Vector{Int}(undef, length(slack_indices))
    g = Vector{Int}(undef, length(variable_indices))
    #Search gb for reducers with minimal slack and add them to the `steps` list
    for row in gb
        copyto!(slacks, 1, row, first_slack_idx, length(slack_indices))
        copyto!(g, 1, row, variable_indices[1], length(variable_indices))
        if IPGBs.GBElements.reduces(g, x)
            steps = new_step_list(steps, slacks)
        end
    end
    #If no potential reducers of x were found, simply return an empty list
    if steps[1] == fill(typemax(Int64), num_slacks)
        return Vector{Int}[]
    end
    return steps
end

"""
    gbelems_with_positive_slack(
    gb::Vector{Vector{Int}},
    slack_indices::Vector{Int}
)::Vector{Vector{Int}}

Return the vectors in `gb` with no negative value and at least one strictly
positive value in the `slack_indices` variables.
"""
function gbelems_with_positive_slack(
    gb::Vector{Vector{Int}},
    slack_indices::Vector{Int}
)::Vector{Vector{Int}}
    if isempty(gb)
        return gb
    end
    pos_slacks = Vector{Int}[]
    for g in gb
        num_zeros = 0
        is_nonnegative = true
        for j in slack_indices
            if g[j] < 0
                is_nonnegative = false
                break
            elseif g[j] == 0
                num_zeros += 1
            end
        end
        if is_nonnegative && num_zeros < length(slack_indices)
            push!(pos_slacks, g)
        end
    end
    return pos_slacks
end

mutable struct MOIPGBState
    instance :: IPInstance #Original multiobjective instance
    ip :: IPInstance #epsilon-constraint IP model for the current problem
    solver::String
    identifier::String

    objective_index::Int #Index of the current objective in `instance`

    current_solution::Vector{Int}
    stats::MultiObjectiveStats.Stats
    test_set::BinomialSet
    test_set_iter::Int
    efficient_points :: Vector{Vector{Int}}
    nondominated_points :: Set{Vector{Int}}
    ideal :: Vector{Int} #(Lower bound of) the ideal point
    nadir :: Vector{Int} #(Upper bound of) the nadir point

    # Variables are classified between:
    original_variable_indices :: Vector{Int} # those that appear in instance
    epsilon_slack_indices :: Vector{Int} # slacks of epsilon constraints

    function MOIPGBState(
        instance, ip, i, init_sol, solver, ideal, nadir,
        vars, slacks, efficient = Vector{Int}[],
        nondominated = Set{Vector{Int}}(), proj_name = nothing, stats = nothing
    )
        if isnothing(stats)
            stats = Stats()
        end
        if isnothing(proj_name)
            proj_name = randstring(RandomDevice(), 12) * "_" * string(i)
        else
            proj_prefix = split(proj_name, "_")[1]
            proj_name = proj_prefix * "_" * string(i)
        end
        new(instance, ip, solver, proj_name, i, copy(init_sol), stats,
            BinomialSet(Vector{Int}[], ip), 0, efficient,
            nondominated, ideal, nadir, vars, slacks
        )
    end
end

#Getters for the most useful data in the MOIPGBState
constraint_matrix(state::MOIPGBState) = state.ip.A
rhs(state::MOIPGBState) = state.ip.b
objective_function(state::MOIPGBState, objective_index::Int) =
    state.ip.C[objective_index, :]
objective_function(state::MOIPGBState) =
    objective_function(state, state.objective_index)
num_variables(state::MOIPGBState) = state.ip.n
num_constraints(state::MOIPGBState) = state.ip.m

function initialize(
    instance :: IPInstance,
    initial_solution :: Vector{Int},
    solver :: String
) :: MOIPGBState
    ideal, nadir = ideal_and_nadir(instance, STD_CLASSICAL_SOLVER)
    vars = collect(1:instance.n)
    slacks = Int[]
    return MOIPGBState(
        instance, copy(instance), 1, initial_solution, solver, ideal, nadir, vars,
        slacks
    )
end

function terminate(state :: MOIPGBState)
    MultiObjectiveStats.terminate(state.stats, state.nondominated_points)
    return state.efficient_points, state.nondominated_points, state.stats
end

function generate_epsilon_constraints(
    model :: JuMP.Model,
    x :: Vector{JuMP.VariableRef},
    prev_state :: MOIPGBState,
    objective_index :: Int
)
    #Add epsilon constraints to all objectives of lower index than
    #objective_index
    epsilon_vars = Int[]
    for i in 1:objective_index-1
        @constraint(model, prev_state.instance.C[i, :]' * x <= prev_state.nadir[i])
        push!(epsilon_vars, length(x) + i)
    end
    return epsilon_vars
end

function generate_ideal_bounds(
    model :: JuMP.Model,
    x :: Vector{JuMP.VariableRef},
    prev_state :: MOIPGBState,
    objective_index :: Int
)
    for i in 1:objective_index-1
        @constraint(model, prev_state.instance.C[i, :]' * x >= prev_state.ideal[i])
    end
end

function generate_nadir_bounds(
    model :: JuMP.Model,
    x :: Vector{JuMP.VariableRef},
    prev_state :: MOIPGBState,
    objective_index :: Int
)
    #Add nadir bound for the current objective for efficiency.
    #GB truncation can eliminate some useless S-vectors when an objective
    #bound given by the nadir point is known.
    i = objective_index
    @constraint(model, prev_state.instance.C[i, :]' * x <= prev_state.nadir[i])
end

function next_objective(
    state :: MOIPGBState,
    initial_solution :: Vector{Int},
    new_objective :: Int,
    use_nadir_bounds :: Bool = true
) :: MOIPGBState
    # Take the initial ip from the instance again
    # Then order its objectives accordingly and add relevant epsilon constraints
    # Finally, add ideal and nadir bounds as needed
    model, x = IPInstances.jump_model(state.instance, optimizer=nothing)
    C = IPInstances.new_first_objective(state.instance, new_objective)
    @objective(model, Min, C * x)
    orig_vars = collect(1:length(x))
    epsilon_vars = generate_epsilon_constraints(model, x, state, new_objective)
    #generate_ideal_bounds(model, x, state, new_objective)
    if use_nadir_bounds
        generate_nadir_bounds(model, x, state, new_objective)
    end
    ip = IPInstance(model)
    initial_solution = IPInstances.extend_feasible_solution(ip, initial_solution)
    return MOIPGBState(
        state.instance, ip, new_objective, initial_solution, state.solver,
        state.ideal, state.nadir, orig_vars, epsilon_vars,
        state.efficient_points, state.nondominated_points, state.identifier,
        state.stats
    )
end

function epsilon_rhs(
    state :: MOIPGBState,
    efficient_point :: Vector{Int},
    step :: Vector{Int}
) :: Vector{Int}
    #Only the first n columns are relevant for the epsilon constraints
    n = length(efficient_point)
    #The last p rows of the constraint matrix are epsilon constraints
    original_num_rows = size(state.instance.A, 1)
    p = length(step)
    #TODO: This should probably be made more flexible later
    epsilon_rows = (original_num_rows+1):(original_num_rows+p)
    C = constraint_matrix(state)[epsilon_rows, 1:n]
    epsilons = C * efficient_point + step
    #TODO: Check this carefully
    extra_slacks = state.ip.b[(original_num_rows+p+1):end]
    return vcat(rhs(state)[1:original_num_rows], epsilons, extra_slacks)
end

function set_current_solution!(
    state :: MOIPGBState,
    efficient_point :: Vector{Int},
    step :: Vector{Int}
)
    epsilons = epsilon_rhs(state, efficient_point, step)
    full_sol = MatrixTools.lift_partial_solution(efficient_point, epsilons, state.ip.A)
    state.current_solution = full_sol
end

"""
    push_efficient_point!(state :: MOIPGBState, point :: Vector{Int})

Update `state`'s efficient and nondominated sets with `point` whenever
`point` is a previously unknown efficient solution.
"""
function push_efficient_point!(state :: MOIPGBState, point :: Vector{Int})
    nondominated = state.instance.C * point
    if !(nondominated in state.nondominated_points)
        push!(state.efficient_points, point)
        push!(state.nondominated_points, nondominated)
    end
end

const TEST_SET_SOLVERS = ["4ti2", "IPGBs"]
const CLASSICAL_SOLVERS = ["GLPK", "CPLEX", "Gurobi", "Cbc"]
const STD_CLASSICAL_SOLVER = IPGBs.DEFAULT_SOLVER
use_test_sets(s::MOIPGBState) = s.solver in TEST_SET_SOLVERS
updated_test_set(s::MOIPGBState) = s.test_set_iter == s.objective_index
has_test_set(s::MOIPGBState) = !isempty(s.test_set) && updated_test_set(s)
compute_new_solution(s :: MOIPGBState, had_test_set :: Bool) =
    s.objective_index == 1 || had_test_set

"""
    moip_gb(s :: MOIPGBState) :: BinomialSet{Vector{Int}, MonomialOrder}

Gröbner Basis / test set for the multi-objective problem corresponding to the
state `s`.

The GB can be computed using either 4ti2 or IPGBs, depending on the solver
specified in `s`.
"""
function moip_gb(
    s :: MOIPGBState;
    gb_size_limit :: Int = 0
) :: BinomialSet{Vector{Int},MonomialOrder}
    @assert use_test_sets(s)
    if s.solver == "4ti2"
        gb = binomials(groebner(s.ip, project_name=s.identifier))
    else #s.solver == "IPGBs"
        gb = groebner_basis(s.ip, gb_size_limit=gb_size_limit)
    end
    @debug "Computed new Gröbner Basis" length(gb)
    return BinomialSet(gb, s.ip)
end

"""
    solve_current_moip!(s::MOIPGBState)

Solve the scalarization of a MOIP given by `s`. If `s.solver` uses test sets,
then a test set is computed, in case one is not already available. Otherwise,
`s.solver` is used directly to solve the scalarization.
"""
function solve_current_moip!(
    s :: MOIPGBState;
    gb_size_limit :: Int = 0
)
    had_test_set = !use_test_sets(s) || has_test_set(s)
    if !use_test_sets(s) #Use a conventional MIP solver
        x, ts, _, _, _ = @timed IPInstances.solve(s.ip)
        efficient_candidate = x[s.original_variable_indices]
        s.stats.solver_time += ts
        s.stats.num_ips += 1
        push_efficient_point!(s, efficient_candidate)
        return x
    elseif !has_test_set(s) #&& use_test_sets(s)
        #Compute a test set and store it in s.test_set
        test_set, time_test_set, _, _, _ = @timed moip_gb(s, gb_size_limit=gb_size_limit)
        s.stats.gb_total_time += time_test_set
        new_gb_size_and_time(s.stats, length(test_set), time_test_set)
        s.test_set = test_set
        s.test_set_iter = s.objective_index
    end
    if compute_new_solution(s, had_test_set)
        #Now use the test set in s to efficiently solve the problem
        _, tnf, _, _, _ = @timed BinomialSets.reduce!(
            s.current_solution, s.test_set, is_monomial_reduction=true
        )
        s.stats.normalform_time += tnf
        s.stats.num_ips += 1
        #Remove slack variables corresponding to the epsilon constraints and
        #push the new efficient solution without them to the pareto set
        efficient_candidate = s.current_solution[s.original_variable_indices]
        push_efficient_point!(s, efficient_candidate)
    end
end

"""
    moip_gb_solve(
    instance::IPInstance,
    initial_solution::Vector{Int};
    solver::String="4ti2"
)::Tuple{Vector{Vector{Int}}, Set{Vector{Int}}, Stats}

Pareto front for multi-objective `instance` given `initial_solution`, using `solver`.

If `solver` is a test-set based solver, applies the project-and-lift algorithm to
obtain a solution.
"""
function moip_gb_solve(
    instance::IPInstance,
    initial_solution::Vector{Int};
    solver::String="IPGBs",
    use_nadir_bounds::Bool=false,
    gb_size_limit::Int=0
)::Tuple{Vector{Vector{Int}}, Set{Vector{Int}}, Stats}
    @debug "Starting multiobjective GB algorithm"
    state = initialize(instance, initial_solution, solver)
    for objective in 1:num_objectives(instance)
        @debug "Starting computation for objective $objective"
        state = next_objective(state, initial_solution, objective, use_nadir_bounds)
        solve_current_moip!(state, gb_size_limit=gb_size_limit)
        pos_slacks = gbelems_with_positive_slack(binomials(state.test_set), state.epsilon_slack_indices)
        for efficient_point in state.efficient_points
            steps = minimal_steps(pos_slacks, efficient_point, state.original_variable_indices, state.epsilon_slack_indices)
            for step in steps
                set_current_solution!(state, efficient_point, step)
                solve_current_moip!(state)
            end
        end
    end
    return terminate(state)
end

function moip_gb_solve(
    instance :: IPInstance;
    solver :: String = "IPGBs",
    use_nadir_bounds :: Bool = false,
    gb_size_limit :: Int = 0
)
    initial_solution = IPInstances.initial_solution(instance)
    return moip_gb_solve(
        instance,
        initial_solution,
        solver=solver,
        use_nadir_bounds=use_nadir_bounds,
        gb_size_limit=gb_size_limit
    )
end

function moip_gb_solve(
    filepath :: String;
    solver :: String = "IPGBs",
    use_nadir_bounds :: Bool = false,
    gb_size_limit :: Int = 0
)
    instance = multiobjective_from_file(filepath)
    return moip_gb_solve(
        instance,
        solver=solver,
        use_nadir_bounds=use_nadir_bounds,
        gb_size_limit=gb_size_limit
    )
end

function moip_walkback(
    instance :: IPInstance
)::Tuple{Vector{Vector{Int}}, Set{Vector{Int}}, Stats}
    #Enumerate all feasible solutions for this instance
    feasibles = enumerate_solutions(instance)
    #Then keep adding them into a nondominated set and checking if it remains
    #nondominated. If it doesn't, remove the dominated solutions from the set
    feasible_vector = [ x for x in feasibles ]
    nondominated_points = Vector{Int}[]
    efficient_points = Vector{Int}[]
    for x in feasible_vector
        y = round.(Int, instance.C * x)
        if is_nondominated(y, nondominated_points, maximization=false)
            remove_dominated!(nondominated_points, efficient_points, y, maximization=false)
            push!(nondominated_points, y)
            push!(efficient_points, x)
        end
    end
    pareto = Set{Vector{Int}}()
    for nondominated in nondominated_points
        push!(pareto, nondominated)
    end
    return efficient_points, pareto, Stats()
end

function moip_walkback(filepath :: String)
    instance = multiobjective_from_file(filepath)
    return moip_walkback(instance)
end

end
