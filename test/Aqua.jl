using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(IPGBs, ambiguities=false)
    Aqua.test_ambiguities(IPGBs)
end
