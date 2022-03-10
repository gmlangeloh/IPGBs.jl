using IPGBs
using IPGBs.FGLM
using IPGBs.FourTi2
using IPGBs.IPInstances
using IPGBs.Markov
using IPGBs.Orders

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
            #TODO I can generalize this into a function that receives an
            #instance later
            instance = IPInstances.random_ipinstance(m, n)
            group_instance = group_relaxation(instance)
            lex_gb = lex_groebner_basis(group_instance)
            det = lex_determinant(lex_gb)
            target_order = MonomialOrder(group_instance,
                                         group_instance.nonnegative_end)
            fglm_gb, t_fglm, _, _, _ = @timed fglm(lex_gb, target_order)
            fourti2_gb, t_4ti2, _, _, _ = @timed groebner(group_instance)
            print(n, " ", m, " ", rep, " random ", det, " ")
            print(t_fglm, " ", t_4ti2, " ")
            println(length(fglm_gb), " ", size(fourti2_gb, 1))
            flush(stdout)
        end
    end
end
