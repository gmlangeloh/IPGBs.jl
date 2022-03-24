using IPGBs
using IPGBs.FGLM
using IPGBs.FourTi2
using IPGBs.IPInstances
using IPGBs.Markov
using IPGBs.Orders

using Random

Random.seed!(0)

instance = IPInstances.random_ipinstance(3, 5)
group_instance = group_relaxation(instance)
lex_gb = lex_groebner_basis(group_instance)
target_order = MonomialOrder(group_instance, group_instance.nonnegative_end)
final_gb = fglm(lex_gb, target_order)
@show final_gb
@show [ lift_vector(v, group_instance) for v in final_gb.basis ]

#Testing the comparison to 4ti2
gb = groebner(group_instance)
@show gb
