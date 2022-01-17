# Unit tests for api/md_result.jl

struct MockResult1 <: MolecularDynamicsResult end
struct MockResult2 <: MolecularDynamicsResult end
Atomistic.kinetic_energy(::MockResult2, t::Integer = 0) = t * u"bohr"
Atomistic.potential_energy(::MockResult2, t::Integer = 0) = -2t * u"bohr"

@testset "api/md_result.jl" begin
    target1 = MockResult1()
    target2 = MockResult2()

    @test_throws Atomistic.UnimplementedError{MockResult1} get_system(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} get_system(target1, 1)
    @test_throws Atomistic.UnimplementedError{MockResult1} get_time_range(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} temperature(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} temperature(target1, 1)
    @test ismissing(reference_temperature(target1))
    @test_throws Atomistic.UnimplementedError{MockResult1} kinetic_energy(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} kinetic_energy(target1, 1)
    @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.potential_energy(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.potential_energy(target1, 1)
    @test total_energy(target2) == 0u"bohr"
    @test total_energy(target2, 1) == -1u"bohr"
    @test_throws Atomistic.UnimplementedError{MockResult1} rdf(target1)
    @test_throws Atomistic.UnimplementedError{MockResult1} rdf(target1, 0.5)
end
