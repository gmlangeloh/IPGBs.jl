"""
Reproduces experiments from Jimenez-Tafur's thesis.
"""

import Random
using Logging

using Knapsack
using MultiObjectiveAlgorithms

function printheader()
    println("family seed n m bin pareto t12 t21 ips totaltime timet12 timet21 timegb timenf timesolver")
end

function runsolve(
    family :: String,
    n :: Int,
    p :: Int,
    repetitions :: Int,
    binary :: Bool
)
    m = 1
    if family == "A"
        generator = knapsack_A
    elseif family == "C"
        generator = knapsack_C
    else
        error("Family has to be either A or C")
    end
    for seed in 1:repetitions
        Random.seed!(seed)
        knapsack = generator(n, objectives=p, binary=binary, format="4ti2")
        kinit = knapsack_initial(knapsack)
        pareto, stats = solve(knapsack, kinit)
        print(family, " ", seed, " ", n, " ", m, " ", binary, " ")
        println(stats)
        flush(stdout)
    end
end

logfile = open("knapsack.log", "w")
logger = SimpleLogger(logfile)
global_logger(logger)

if length(ARGS) < 1
    error("Missing command line argument")
end

table = ARGS[1]

if table == "4.15"
    #Table 4.15: 3 objectives, binary
    printheader()
    runsolve("A", 10, 3, 30, true)
    runsolve("A", 15, 3, 30, true)
elseif table == "4.16"
    #Table 4.16: 3 objectives, unbounded
    printheader()
    runsolve("A", 25, 3, 30, false)
    runsolve("A", 50, 3, 30, false)
    runsolve("A", 75, 3, 30, false)
    runsolve("A", 100, 3, 30, false)
    runsolve("C", 25, 3, 30, false)
    runsolve("C", 50, 3, 30, false)
    runsolve("C", 75, 3, 30, false)
    runsolve("C", 100, 3, 30, false)
elseif table == "4.17"
    #Table 4.17: 4 objectives, unbounded
    printheader()
    runsolve("A", 25, 4, 30, false)
    runsolve("A", 50, 4, 30, false)
    runsolve("A", 75, 4, 30, false)
    runsolve("C", 25, 4, 30, false)
    runsolve("C", 50, 4, 30, false)
    runsolve("C", 75, 4, 30, false)
elseif table == "4.18"
    #Table 4.18: 5 objectives, unbounded
    printheader()
    runsolve("A", 25, 5, 30, false)
    runsolve("A", 50, 5, 30, false)
    runsolve("A", 75, 5, 30, false)
    runsolve("C", 25, 5, 30, false)
    runsolve("C", 50, 5, 30, false)
    runsolve("C", 75, 5, 30, false)
else
    error("Unknown argument")
end
