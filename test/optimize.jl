using CPLEX
using IPGBs
using JuMP
using MIPMatrixTools
using MIPMatrixTools.IPInstances

function test_optimize(path; init_sol = nothing)
    model = read_from_file(path)
    set_optimizer(model, CPLEX.Optimizer)
    set_silent(model)
    optimize!(model)
    value = round(Int, objective_value(model))
    if isnothing(init_sol)
        ret, t, _, _, _ = @timed IPGBs.optimize(path, quiet=false)
    else
        ret, t, _, _, _ = @timed IPGBs.optimize(path, quiet=false, solution=init_sol)
    end
    val = ret[2]
    println("IPGBs opt: ", val, " time: ", t)
    println("Cplex opt: ", value)
    return ret[1]
end
