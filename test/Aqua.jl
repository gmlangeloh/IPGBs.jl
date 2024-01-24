using Aqua

@testset "Aqua.jl" begin
    println("Aqua.jl tests")
    Aqua.test_all(IPGBs, ambiguities=false)
    Aqua.test_ambiguities(IPGBs)
end
