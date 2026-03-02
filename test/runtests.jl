using OrbitPropagationLibrarySOFA
using OrbitPropagationLibraryNRLMSIS
using Test

const tests = [
    "nasa",
]
@testset "OrbitPropagationLibraryNRLMSIS.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
