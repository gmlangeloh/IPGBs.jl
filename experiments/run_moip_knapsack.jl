"""
Runs an epsilon-constraint based bi-objective solver over randomly generated
bi-objective knapsack instances.

Parameters:
- a string containing the desired families of instances (a substring of "ABCD")
or tables to reproduce (a substring of "123456")
- number of items for the knapsacks
- number of random instances of each family to be generated

If the number of a table is given, no other parameters are necessary.
"""

using MultiObjectiveInstances
using MultiObjectiveInstances.Knapsack

using IPGBs
using IPGBs.MultiObjectiveStats
using IPGBs.MultiObjectiveTools

using Logging
using Random

using CPLEX
using Cbc

function printheader(algorithm, num_objectives :: Int)
    print("family seed n m bin ")
    if algorithm == "grobnecon"
        MultiObjectiveStats.header(num_objectives)
    elseif algorithm == "augmecon"
        #TODO: Use the same header in augmecon and a general solver
        println("pareto totaltime timepayoff timesolver timemodel timenadir")
    else
        println("pareto totaltime solvertime nummodels")
    end
end

function runsolve(
    family :: String,
    n :: Int,
    m :: Int,
    p :: Int,
    repetitions :: Int;
    binary :: Bool = false,
    solver :: String = "4ti2",
    algorithm :: String = "grobnecon"
)
    @assert family in ["A", "B", "C", "D"]
    printheader(algorithm, p)
    generator = knapsack_A
    if family == "A"
        generator = knapsack_A
    elseif family == "B"
        generator = knapsack_B
    elseif family == "C"
        generator = knapsack_C
    elseif family == "D"
        generator = knapsack_D
    end
    for seed in 1:repetitions
        Random.seed!(seed)
        knapsack = generator(n, m, objectives = p, binary = binary)
        if algorithm == "grobnecon"
            knapsack = fourti2_stdform(knapsack)
            filename = "knapsack_$(n)_$(m)_$(p)_$(seed).dat"
            path = "./moip_instances/"
            fullname = path * filename
            MultiObjectiveInstances.write_to_file(knapsack, fullname)
            kinit = knapsack_initial(knapsack)
            _, _, stats = moip_gb_solve(knapsack, kinit, solver = solver)
            print(family, " ", seed, " ", n, " ", m, " ", binary, " ")
            println(stats)
        end
        flush(stdout)
    end
end

if length(ARGS) > 0
    if length(ARGS) < 2
        error("Missing command line arguments.")
    end

    algorithm = ARGS[1]
    if !(algorithm in [ "grobnecon" ])
        error("Algorithm must be grobnecon.")
    end
    families = ARGS[2]
    #Set up the log file for this run
    logfile = open("knapsack.log", "w")
    logger = SimpleLogger(logfile, Logging.Debug)
    global_logger(logger)

    #Is the current code trying to reproduce some tables from Hartillo-Hermoso
    #et al (2019)?
    finished = false
    if families == "1"
        runsolve("A", 50, 1, 2, 30, algorithm=algorithm)
        runsolve("B", 50, 1, 2, 30, algorithm=algorithm)
        runsolve("C", 50, 1, 2, 30, algorithm=algorithm)
        runsolve("D", 50, 1, 2, 30, algorithm=algorithm)
        finished = true
    elseif families == "2"
        runsolve("D", 60, 1, 2, 30, algorithm=algorithm)
        runsolve("D", 70, 1, 2, 30, algorithm=algorithm)
        runsolve("D", 80, 1, 2, 30, algorithm=algorithm)
        runsolve("D", 90, 1, 2, 30, algorithm=algorithm)
        runsolve("D", 100, 1, 2, 30, algorithm=algorithm)
        finished = true
    elseif families == "3"
        runsolve("A", 10, 1, 2, 30, binary = true, algorithm=algorithm)
        runsolve("A", 15, 1, 2, 30, binary = true, algorithm=algorithm)
        runsolve("A", 20, 1, 2, 30, binary = true, algorithm=algorithm)
        runsolve("A", 25, 1, 2, 30, binary = true, algorithm=algorithm)
        runsolve("A", 30, 1, 2, 30, binary = true, algorithm=algorithm)
        finished = true
    elseif families == "4"
        runsolve("A", 20, 2, 2, 30, algorithm=algorithm)
        runsolve("A", 25, 2, 2, 30, algorithm=algorithm)
        runsolve("A", 30, 2, 2, 30, algorithm=algorithm)
        runsolve("A", 35, 2, 2, 30, algorithm=algorithm)
        runsolve("A", 40, 2, 2, 30, algorithm=algorithm)
        finished = true
    elseif families == "5"
        runsolve("A", 20, 3, 2, 30, algorithm=algorithm)
        runsolve("A", 25, 3, 2, 30, algorithm=algorithm)
        runsolve("A", 30, 3, 2, 30, algorithm=algorithm)
        runsolve("A", 35, 3, 2, 30, algorithm=algorithm)
        finished = true
    elseif families == "6"
        runsolve("A", 100, 1, 2, 30, algorithm=algorithm)
        runsolve("A", 200, 1, 2, 30, algorithm=algorithm)
        runsolve("A", 300, 1, 2, 30, algorithm=algorithm)
        runsolve("A", 400, 1, 2, 30, algorithm=algorithm)
        runsolve("A", 500, 1, 2, 30, algorithm=algorithm)
        finished = true
    end

    if finished
        exit()
    end

    if length(ARGS) < 3
        error("Missing command line arguments.")
    end

    n = parse(Int, ARGS[3])

    if length(ARGS) >= 4
        objectives = parse(Int, ARGS[4])
    else
        objectives = 2
    end
    if length(ARGS) >= 5
        repetitions = parse(Int, ARGS[5])
    else
        repetitions = 1
    end
    if length(ARGS) >= 6
        solver = ARGS[6]
    else
        solver = "4ti2"
    end
    binary = false
    if length(ARGS) >= 7
        if ARGS[7] == "binary"
            binary = true
        end
    end
    if occursin("A", families)
        runsolve("A", n, 1, objectives, repetitions, solver=solver, binary=binary,
                 algorithm=algorithm)
    end
    if occursin("B", families)
        runsolve("B", n, 1, objectives, repetitions, solver=solver, binary=binary,
                 algorithm=algorithm)
    end
    if occursin("C", families)
        runsolve("C", n, 1, objectives, repetitions, solver=solver, binary=binary,
                 algorithm=algorithm)
    end
    if occursin("D", families)
        runsolve("D", n, 1, objectives, repetitions, solver=solver, binary=binary,
                 algorithm=algorithm)
    end
end
