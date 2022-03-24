using IPGBs
using IPGBs.FGLM
using IPGBs.FourTi2
using IPGBs.IPInstances
using IPGBs.Markov
using IPGBs.Orders

using MultiObjectiveInstances
using Random

#TODO this should be somewhere else, but I'll leave it here for now. The
#algorithms themselves don't need this yet
function lex_determinant(lex_gb)
    n = length(lex_gb[1])
    det = 1
    @assert length(lex_gb) == n
    for i in 1:n
        g = lex_gb[i]
        det *= g[i]
    end
    return det
end

function run_instance(instance)
    group_instance = group_relaxation(instance)
    println(group_instance.model)
    lex_gb = lex_groebner_basis(group_instance)
    det = lex_determinant(lex_gb)
    target_order = MonomialOrder(group_instance,
        group_instance.nonnegative_end)
    fglm_gb, t_fglm, _, _, _ = @timed fglm(lex_gb, target_order)
    fourti2_gb, t_4ti2, _, _, _ = @timed groebner(group_instance)
    print(det, " ", t_fglm, " ", t_4ti2, " ")
    println(length(fglm_gb), " ", size(fourti2_gb, 1))
    flush(stdout)
end

function run_random_instances()
    Random.seed!(0)
    ns = collect(5:25)
    ms = collect(1:20)
    println("n m rep type det time_fglm time_4ti2 len_fglm len_4ti2")
    for n in ns
        for m in ms
            if n <= m
                continue
            end
            for rep in 1:10
                instance = IPInstances.random_ipinstance(m, n)
                print(n, " ", m, " ", rep, " random ")
                run_instance(instance)
            end
        end
    end
end

function run_knapsacks()
    Random.seed!(0)
    ns = collect(5:20)
    for n in ns
        for rep in 1:10
            knapsack = MultiObjectiveInstances.Knapsack.knapsack_A(
                n, format = "4ti2", binary = true)
            u = ones(Int, size(knapsack.A, 2))
            instance = IPInstances.IPInstance(
                knapsack.A, knapsack.b, knapsack.C, u,
                apply_normalization = false, invert_objective = false
            )
            print(n, " ", n + 1, " ", rep, " knapsack ")
            run_instance(instance)
        end
    end
end

run_knapsacks()
