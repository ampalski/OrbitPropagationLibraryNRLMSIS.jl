using OrbitPropagationLibrarySOFA
using OrbitPropagationLibraryNRLMSIS
using Test
using Aqua

const tests = [
    "nasa",
    "aqua",
]
@testset "OrbitPropagationLibraryNRLMSIS.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
