# Unit tests for implementations/molly/molly_base.jl

@testset "implementations/molly/molly_base.jl" begin
    particles = [
        AtomsBase.Atom(:Ar, [-2.0, -1.0, 3.0]u"bohr", [3.0, 5.0, 21.0]u"bohr * hartree / ħ_au"; a=:b),
        AtomsBase.Atom(:Ar, [-7.0, 2.0, 13.0]u"bohr", [-3.0, 7.0, 1.0]u"bohr * hartree / ħ_au"; c=:d)
    ]
    box = (@SVector [(@SVector [1.5, 0.0, 0.0]), (@SVector [0.0, 1.5, 0.0]), (@SVector [0.0, 0.0, 1.5])])u"bohr"
    boundary_conditions = @SVector [Periodic(), Periodic(), Periodic()]
    flexible = FlexibleSystem(particles, box, boundary_conditions)
    fast = FastSystem(particles, box, boundary_conditions)

    system1 = System(flexible)
    system2 = System(fast)
    system3 = System(system1; loggers=Dict("c" => CoordinateLogger(1)))

    @test system1.atoms == [
        Molly.Atom(; index=1, mass=auconvert(39.9481u"u"))
        Molly.Atom(; index=2, mass=auconvert(39.9481u"u"))
    ]
    @test system1.atoms_data == [
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}(:a => :b)),
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}(:c => :d))
    ]
    @test system1.coords == [
        (@SVector [-2.0, -1.0, 3.0])u"bohr",
        (@SVector [-7.0, 2.0, 13.0])u"bohr"
    ]
    @test system1.velocities == [
        (@SVector [3.0, 5.0, 21.0])u"bohr * hartree / ħ_au",
        (@SVector [-3.0, 7.0, 1.0])u"bohr * hartree / ħ_au"
    ]
    @test system1.box_size == (@SVector [1.5, 1.5, 1.5])u"bohr"
    @test system1.force_units == u"hartree / bohr"
    @test system1.energy_units == u"hartree"

    @test system2.atoms == [
        Molly.Atom(; index=1, mass=auconvert(39.9481u"u"))
        Molly.Atom(; index=2, mass=auconvert(39.9481u"u"))
    ]
    @test system2.atoms_data == [
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}()),
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}())
    ]
    @test system2.coords == [
        (@SVector [-2.0, -1.0, 3.0])u"bohr",
        (@SVector [-7.0, 2.0, 13.0])u"bohr"
    ]
    @test system2.velocities == [
        (@SVector [0.0, 0.0, 0.0])u"bohr * hartree / ħ_au",
        (@SVector [0.0, 0.0, 0.0])u"bohr * hartree / ħ_au"
    ]
    @test system2.box_size == (@SVector [1.5, 1.5, 1.5])u"bohr"
    @test system2.force_units == u"hartree / bohr"
    @test system2.energy_units == u"hartree"

    @test system3.atoms == [
        Molly.Atom(; index=1, mass=auconvert(39.9481u"u"))
        Molly.Atom(; index=2, mass=auconvert(39.9481u"u"))
    ]
    @test system3.atoms_data == [
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}(:a => :b)),
        Atomistic.AugmentedAtomData(:Ar, Dict{Symbol,Any}(:c => :d))
    ]
    @test system3.coords == [
        (@SVector [-2.0, -1.0, 3.0])u"bohr",
        (@SVector [-7.0, 2.0, 13.0])u"bohr"
    ]
    @test system3.velocities == [
        (@SVector [3.0, 5.0, 21.0])u"bohr * hartree / ħ_au",
        (@SVector [-3.0, 7.0, 1.0])u"bohr * hartree / ħ_au"
    ]
    @test system3.box_size == (@SVector [1.5, 1.5, 1.5])u"bohr"
    @test system3.force_units == u"hartree / bohr"
    @test system3.energy_units == u"hartree"
    @test system3.loggers["c"] isa typeof(CoordinateLogger(1))

    data = Atomistic.AugmentedAtomData(:Ar, a=:b)

    @test hasproperty(data, :element)
    @test hasproperty(data, :a)
    @test !hasproperty(data, :b)
    @test data.element == :Ar
    @test data.a == :b
    @test_throws KeyError data.b
    @test propertynames(data) == (:element, :a)
    @test propertynames(data, true) == (:element, :data, :a)

    system = System(;
        atoms=[Molly.Atom(; index=1, mass=auconvert(39.9481u"u"))],
        atoms_data=[data],
        coords=[(@SVector [-2.0, -1.0, 3.0])u"bohr"],
        box_size=1.5u"bohr"
    )
    view = system[1]

    @test hasproperty(view, :a)
    @test !hasproperty(view, :b)
    @test view.a == :b
    @test_throws KeyError view.b
    @test propertynames(view) == (:system, :index, :a)
    @test propertynames(view, true) == (:system, :index, :a)

    atom = AtomsBase.Atom(data, (@SVector [-2.0, -1.0, 3.0])u"bohr", (@SVector [3.0, 5.0, 21.0])u"bohr * hartree / ħ_au")

    @test atom == particles[1]

    @test AtomsBase.atomic_symbol(system1, 1) == :Ar
end
