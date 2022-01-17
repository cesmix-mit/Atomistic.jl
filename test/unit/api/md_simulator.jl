# Unit tests for api/md_simulator.jl

struct MockSystem <: AbstractSystem{3} end
struct MockSimulator <: MolecularDynamicsSimulator end
struct MockPotential <: InteratomicPotentials.ArbitraryPotential end

@testset "api/md_simulator.jl" begin
    @test_throws Atomistic.UnimplementedError{MockSimulator} simulate(MockSystem(), MockSimulator(), MockPotential())
end
