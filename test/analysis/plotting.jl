# Unit tests for analysis/plotting.jl

@testset "analysis/plotting.jl" begin
    struct MockTimeResult <: MolecularDynamicsResult end
    Atomistic.get_time_range(::MockTimeResult) = [0, 1, 2, 4, 5, 7, 8]u"ns"

    target = MockTimeResult()

    @test Atomistic.plotting_time_range(target, 1) == [(1, 0u"ps"), (2, 1000u"ps"), (3, 2000u"ps"), (4, 4000u"ps"), (5, 5000u"ps"), (6, 7000u"ps"), (7, 8000u"ps")]
    @test Atomistic.plotting_time_range(target, 2) == [(1, 0u"ps"), (3, 2000u"ps"), (5, 5000u"ps"), (7, 8000u"ps")]
    @test Atomistic.plotting_time_range(target, 6) == [(1, 0u"ps"), (7, 8000u"ps")]
    @test Atomistic.plotting_time_range(target, 7) == [(1, 0u"ps")]
    @test Atomistic.plotting_time_range(target, 8) == [(1, 0u"ps")]
end
