# Unit tests for integrations/atomsbase_integration.jl

@testset "integrations/atomsbase_integration.jl" begin
    struct MockSystem1D <: AbstractSystem{1} end
    AtomsBase.bounding_box(s::MockSystem1D) = infinite_box(3)
    AtomsBase.boundary_conditions(s::MockSystem1D) = [Periodic(), Periodic(), Periodic()]
    Base.length(s::MockSystem1D) = 3
    AtomsBase.species_type(s::MockSystem1D) = Symbol
    Base.getindex(s::MockSystem1D, index::Int) = Symbol(index)
    AtomsBase.position(s::MockSystem1D) = [11, 12, 13]u"m"
    AtomsBase.atomic_symbol(s::MockSystem1D) = [Symbol(11), Symbol(12), Symbol(13)]
    AtomsBase.atomic_number(s::MockSystem1D) = [10, 10, 10]
    AtomsBase.atomic_mass(s::MockSystem1D) = [10, 10, 10]u"kg"
    AtomsBase.velocity(s::MockSystem1D) = missing
    AtomsBase.position(s::MockSystem1D, i) = i * u"m"
    AtomsBase.atomic_symbol(s::MockSystem1D, i) = Symbol(i)
    AtomsBase.atomic_number(s::MockSystem1D, i) = 0
    AtomsBase.atomic_mass(s::MockSystem1D, i) = 0u"kg"
    AtomsBase.velocity(s::MockSystem1D, i) = missing

    target = DynamicSystem(MockSystem1D(), 3.14u"s")

    @test bounding_box(target) == infinite_box(3)
    @test boundary_conditions(target) == [Periodic(), Periodic(), Periodic()]
    @test periodicity(target) == [true, true, true]
    @test n_dimensions(target) == 1
    @test length(target) == 3
    @test size(target) == (3,)
    @test species_type(target) == Symbol
    @test [x for x ∈ target] == [Symbol(1), Symbol(2), Symbol(3)]
    @test position(target) == [11, 12, 13]u"m"
    @test ismissing(velocity(target))
    @test AtomsBase.atomic_symbol(target) == [Symbol(11), Symbol(12), Symbol(13)]
    @test atomic_number(target) == [10, 10, 10]
    @test atomic_mass(target) == [10, 10, 10]u"kg"
    @test position(target, 1) == 1u"m"
    @test ismissing(velocity(target, 1))
    @test AtomsBase.atomic_symbol(target, 1) == Symbol(1)
    @test atomic_number(target, 1) == 0
    @test atomic_mass(target, 1) == 0u"kg"
end
