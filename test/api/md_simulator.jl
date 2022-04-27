# Unit tests for api/md_simulator.jl

@testset "api/md_simulator.jl" begin
    struct MockSystem <: AbstractSystem{3} end
    struct MockSimulator <: MolecularDynamicsSimulator end
    struct MockPotential <: InteratomicPotentials.AbstractPotential end

    @test_throws UnimplementedError{MockSimulator} simulate(MockSystem(), MockSimulator(), MockPotential())
end
