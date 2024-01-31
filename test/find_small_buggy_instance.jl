using MIPMatrixTools.CombinatorialOptimizationInstances
using MIPMatrixTools.IPInstances
using MIPMatrixTools.GBTools
using IPGBs
using IPGBs.FourTi2
using Random
using JuMP

function find_buggy_lap()
    Random.seed!(0)
    for _ in 1:1000
        lap, _ = generate_lap(3)
        ip = IPInstance(lap)
        mygb = groebner_basis(ip)
        theirgb = GBTools.tovector(groebner(ip))
        if !GBTools.isequal(mygb, theirgb)
            println("Found buggy lap")
            write_to_file(lap, "lap_3_2.mps")
            break
        end
    end
end

find_buggy_lap()
