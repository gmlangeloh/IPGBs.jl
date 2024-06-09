#Try generating set packing and cover problems that aren't too easy or hard for GBs

using IPGBs
using IPGBs.FourTi2
using IPGBs.CombinatorialOptimizationInstances
using IPGBs.IPInstances

#Here, m is the number of variables (sets), n is the number of constraints (elements)
function try_packing(n, m, p)
    model = generate_set_packing(n, m, p)
    instance = IPInstance(model)
    println("Generated instance, running IPGBs and 4ti2")
    ipgbs, t_ipgbs, _, _, _ = @timed IPGBs.groebner_basis(instance)
    println("IPGBs time: ", t_ipgbs, " basis size: ", length(ipgbs))
    old_ipgbs, t_old_ipgbs, _, _, _ = @timed IPGBs.groebner_basis(instance, use_binary_truncation=false)
    println("IPGBs (no bin truncation): ", t_old_ipgbs, " basis size: ", length(old_ipgbs))
    fti2, t_fti2, _, _, _ = @timed FourTi2.groebner(instance)
    println("4ti2 time: ", t_fti2, " basis size: ", size(fti2, 1))
end

#Here m is the number of constraints (sets), n is the number of variables (elements)
function try_covering(n, m, p)
    model,_ = generate_set_cover(n, m, p)
    instance = IPInstance(model)
    println("Generated instance, running IPGBs and 4ti2")
    ipgbs, t_ipgbs, _, _, _ = @timed IPGBs.groebner_basis(instance)
    println("IPGBs time: ", t_ipgbs, " basis size: ", length(ipgbs))
    old_ipgbs, t_old_ipgbs, _, _, _ = @timed IPGBs.groebner_basis(instance, use_binary_truncation=false)
    println("IPGBs (no bin truncation): ", t_old_ipgbs, " basis size: ", length(old_ipgbs))
    fti2, t_fti2, _, _, _ = @timed FourTi2.groebner(instance)
    println("4ti2 time: ", t_fti2, " basis size: ", size(fti2, 1))
end
