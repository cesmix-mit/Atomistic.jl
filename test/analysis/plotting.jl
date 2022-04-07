# Unit tests for analysis/plotting.jl

@testset "analysis/plotting.jl" begin
    struct MockPlotResult <: MolecularDynamicsResult end
    Atomistic.get_time_range(::MockPlotResult) = [0, 1, 2]u"ns"

    target = MockPlotResult()

    extractor = Atomistic.extract_property(target, get_time, u"ps")

    @test extractor.(1:length(target)) == [0u"ps", 1000u"ps", 2000u"ps"]
end
