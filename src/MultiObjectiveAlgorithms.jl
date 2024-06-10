"""
A module to solve multi-objective integer programming problems using
epsilon-constraint methods and test sets / Gröbner Bases.
"""
module MultiObjectiveAlgorithms
export moip_gb_solve

using IPGBs
using IPGBs.BinomialSets
using IPGBs.FourTi2
using IPGBs.IPInstances
using IPGBs.MatrixTools
using IPGBs.Orders

using MultiObjectiveInstances

using IPGBs.MultiObjectiveTools
using IPGBs.SingleObjective
using IPGBs.MultiObjectiveStats

using JuMP
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
    instance :: MultiObjectiveInstance #Full multiobjective instance
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
            BinomialSet(Vector{Int}[], ip.C, ip.A, ip.b), 0, efficient,
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

#TODO: I should probably move this function somewhere else
function new_first_objective(
    instance :: MultiObjectiveInstance,
    first_objective :: Int
) :: Array{Int, 2}
    obj_permutation = collect(1:num_objectives(instance))
    obj_permutation[first_objective] = 1
    obj_permutation[1] = first_objective
    return instance.C[obj_permutation, :]
end

"""
    initial_jump_model(
    instance :: MultiObjectiveInstance,
    first_objective :: Int = 1
)

Return a JuMP Model corresponding to the MultiObjectiveInstance.

The MultiObjectiveInstance is assumed to be in the form
min C * x
s.t. Ax = b
x >= 0, x Integer.

The JuMP Model keeps all objectives in its objective function matrix.
"""
function initial_jump_model(
    instance :: MultiObjectiveInstance,
    first_objective :: Int = 1
)
    model = JuMP.Model()
    n = MultiObjectiveInstances.num_variables(instance)
    m = MultiObjectiveInstances.num_constraints(instance)
    @variable(model, x[1:n] >= 0, Int)
    #For now, simply give IPInstances the whole C matrix with all objectives
    C = new_first_objective(instance, first_objective)
    A = instance.A
    b = instance.b
    @objective(model, Min, C * x)
    #Write down the constraints according to the matrix A and rhs b
    for i in 1:m
        @constraint(model, sum(A[i, j] * x[j] for j in 1:n) == b[i])
    end
    #Epsilon constraints will be added later as needed
    return model, x
end

initial_ip_instance(instance, first_objective = 1) =
    IPInstance(initial_jump_model(instance, first_objective)[1])

function initialize(
    instance :: MultiObjectiveInstance,
    initial_solution :: Vector{Int},
    solver :: String
) :: MOIPGBState
    ideal = ideal_point(instance, STD_CLASSICAL_SOLVER)
    nadir = nadir_bound(instance, STD_CLASSICAL_SOLVER)
    ip = initial_ip_instance(instance)
    vars = collect(1:MultiObjectiveInstances.num_variables(instance))
    slacks = Int[]
    return MOIPGBState(
        instance, ip, 1, initial_solution, solver, ideal, nadir, vars, slacks
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
    new_objective :: Int
) :: MOIPGBState
    # Take the initial ip from the instance again
    # Then order its objectives accordingly and add relevant epsilon constraints
    # Finally, add ideal and nadir bounds as needed
    model, x = initial_jump_model(state.instance, new_objective)
    orig_vars = collect(1:length(x))
    epsilon_vars = generate_epsilon_constraints(model, x, state, new_objective)
    #generate_ideal_bounds(model, x, state, new_objective)
    #generate_nadir_bounds(model, x, state, new_objective)
    ip = IPInstance(model)
    return MOIPGBState(
        state.instance, ip, new_objective, initial_solution, state.solver,
        state.ideal, state.nadir, orig_vars, epsilon_vars, state.efficient_points,
        state.nondominated_points, state.identifier, state.stats
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
    p = length(step)
    C = constraint_matrix(state)[end-p+1:end, 1:n]
    epsilons = C * efficient_point + step
    #Return the right-hand side vector with the last p elements set to epsilons
    return vcat(rhs(state)[1:end-p], epsilons)
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
const CLASSICAL_SOLVERS = ["CPLEX", "Gurobi", "Cbc"]
const STD_CLASSICAL_SOLVER = "Cbc"
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
function moip_gb(s::MOIPGBState)::BinomialSet{Vector{Int},MonomialOrder}
    @assert use_test_sets(s)
    if s.solver == "4ti2"
        gb = binomials(groebner(s.ip, project_name=s.identifier))
    else #s.solver == "IPGBs"
        gb = groebner_basis(s.ip)
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
function solve_current_moip!(s::MOIPGBState)
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
        test_set, time_test_set, _, _, _ = @timed moip_gb(s)
        s.stats.gb_total_time += time_test_set
        new_gb_size_and_time(s.stats, length(test_set), time_test_set)
        s.test_set = test_set
        s.test_set_iter = s.objective_index
    end
    if compute_new_solution(s, had_test_set)
        #Now use the test set in s to efficiently solve the problem
        _, tnf, _, _, _ = @timed BinomialSets.reduce!(
            s.current_solution, s.test_set
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
    instance::MultiObjectiveInstance,
    initial_solution::Vector{Int};
    solver::String="4ti2"
)::Tuple{Vector{Vector{Int}},Set{Vector{Int}},Stats}

Pareto front for multi-objective `instance` given `initial_solution`, using `solver`.

If `solver` is a test-set based solver, applies the project-and-lift algorithm to
obtain a solution.
"""
function moip_gb_solve(
    instance::MultiObjectiveInstance,
    initial_solution::Vector{Int};
    solver::String="4ti2"
)::Tuple{Vector{Vector{Int}}, Set{Vector{Int}}, Stats}
    @debug "Starting multiobjective GB algorithm"
    state = initialize(instance, initial_solution, solver)
    for objective in 1:num_objectives(instance)
        @debug "Starting computation for objective $objective"
        state = next_objective(state, initial_solution, objective)
        solve_current_moip!(state)
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

end
