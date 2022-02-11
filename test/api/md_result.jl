# Unit tests for api/md_result.jl

@testset "api/md_result.jl" begin
    @testset "Unimplemented Functions" begin
        struct MockResult1 <: MolecularDynamicsResult end

        target = MockResult1()

        @test_throws Atomistic.UnimplementedError{MockResult1} get_time_range(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_num_bodies(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_bounding_box(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_boundary_conditions(target)

        @test_throws Atomistic.UnimplementedError{MockResult1} get_positions(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_positions(target, 1)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_velocities(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_velocities(target, 1)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_particles(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} get_particles(target, 1)

        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.temperature(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.temperature(target, 1)
        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.kinetic_energy(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.kinetic_energy(target, 1)
        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.potential_energy(target)
        @test_throws Atomistic.UnimplementedError{MockResult1} Atomistic.potential_energy(target, 1)
    end
    @testset "Default Implementations" begin
        struct MockResult2 <: MolecularDynamicsResult end
        Atomistic.get_time_range(::MockResult2) = collect(0:2:100)u"ns"
        Atomistic.get_num_bodies(::MockResult2) = 3
        Atomistic.get_bounding_box(::MockResult2) = (@SVector [(@SVector [100, 0, 0]), (@SVector [0, 100, 0]), (@SVector [0, 0, 100])])u"bohr"
        Atomistic.get_boundary_conditions(::MockResult2) = @SVector [Periodic(), Periodic(), Periodic()]
        Atomistic.get_positions(::MockResult2, t::Integer) = [
            (@SVector [0, 0, 0])u"bohr",
            (@SVector [-t, 0, t])u"bohr",
            (@SVector [t, 0, -t])u"bohr"
        ]
        Atomistic.get_particles(::MockResult2, t::Integer) = [
            AtomsBase.Atom(:Ar, (@SVector [0, 0, 0])u"bohr", (@SVector [0, 0, 0])u"bohr/ns"),
            AtomsBase.Atom(:Ar, (@SVector [-t, 0, t])u"bohr", (@SVector [-2t, 0, 2t])u"bohr/ns"),
            AtomsBase.Atom(:Ar, (@SVector [t, 0, -t])u"bohr", (@SVector [2t, 0, -2t])u"bohr/ns")
        ]
        Atomistic.kinetic_energy(::MockResult2, t::Integer) = t * u"bohr"
        Atomistic.potential_energy(::MockResult2, t::Integer) = -2t * u"bohr"

        atoms_equal(a1::AtomsBase.Atom, a2::AtomsBase.Atom) = propertynames(a1) == propertynames(a2) && all(getproperty(a1, p) == getproperty(a2, p) for p ∈ propertynames(a1))

        target = MockResult2()

        @test ismissing(reference_temperature(target))

        @test get_time(target) == 100u"ns"
        @test get_time(target, 1) == 0u"ns"

        system = get_system(target)
        @test bounding_box(system) == (@SVector [(@SVector [100, 0, 0]), (@SVector [0, 100, 0]), (@SVector [0, 0, 100])])u"bohr"
        @test boundary_conditions(system) == @SVector [Periodic(), Periodic(), Periodic()]
        @test atoms_equal(system[1], AtomsBase.Atom(:Ar, (@SVector [0, 0, 0])u"bohr", (@SVector [0, 0, 0])u"bohr/ns"))
        @test atoms_equal(system[2], AtomsBase.Atom(:Ar, (@SVector [-51, 0, 51])u"bohr", (@SVector [-102, 0, 102])u"bohr/ns"))
        @test atoms_equal(system[3], AtomsBase.Atom(:Ar, (@SVector [51, 0, -51])u"bohr", (@SVector [102, 0, -102])u"bohr/ns"))
        system1 = get_system(target, 1)
        @test bounding_box(system1) == (@SVector [(@SVector [100, 0, 0]), (@SVector [0, 100, 0]), (@SVector [0, 0, 100])])u"bohr"
        @test boundary_conditions(system1) == @SVector [Periodic(), Periodic(), Periodic()]
        @test atoms_equal(system1[1], AtomsBase.Atom(:Ar, (@SVector [0, 0, 0])u"bohr", (@SVector [0, 0, 0])u"bohr/ns"))
        @test atoms_equal(system1[2], AtomsBase.Atom(:Ar, (@SVector [-1, 0, 1])u"bohr", (@SVector [-2, 0, 2])u"bohr/ns"))
        @test atoms_equal(system1[3], AtomsBase.Atom(:Ar, (@SVector [1, 0, -1])u"bohr", (@SVector [2, 0, -2])u"bohr/ns"))

        @test Atomistic.total_energy(target) == -51u"bohr"
        @test Atomistic.total_energy(target, 1) == -1u"bohr"
        rdf1 = Atomistic.rdf(target)
        rdf2 = Atomistic.rdf(target, 1)
        rdf3 = Atomistic.rdf(target, 1, 51)
        @test rdf1 isa NamedTuple{(:r, :g),<:Tuple{<:AbstractVector{<:Unitful.Length},<:AbstractVector{Float64}}}
        @test rdf1 == rdf2
        @test rdf1 == rdf3
    end
end
