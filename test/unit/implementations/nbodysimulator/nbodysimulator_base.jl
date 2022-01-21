# Unit tests for implementations/nbodysimulator/nbodysimulator_base.jl

@testset "implementations/nbodysimulator/nbodysimulator_base.jl" begin
    body1 = Atomistic.ElementMassBody((@SVector [1, 2, 3])u"nm", (@SVector [4, 5, 6])u"nm/s", elements[:Ar]; meta = "data")

    @test body1 isa Atomistic.ElementMassBody{Float64,Float64}
    @test body1.r ≈ (@SVector [1, 2, 3])austrip(1u"nm")
    @test body1.v ≈ (@SVector [4, 5, 6])austrip(1u"nm/s")
    @test body1.m ≈ 39.9481austrip(1u"u")
    @test body1.symbol == :Ar
    @test body1.number == 18
    @test body1.data == Dict(:meta => "data")

    body2 = Atomistic.ElementMassBody(body1, (@SVector [2, 3, 4]), Float32.(@SVector [5, 6, 7]))

    @test body2 isa Atomistic.ElementMassBody{Float32,Float64}
    @test body2.r ≈ @SVector [2, 3, 4]
    @test body2.v ≈ @SVector [5, 6, 7]
    @test body2.m ≈ 39.9481austrip(1u"u")
    @test body2.symbol == :Ar
    @test body2.number == 18
    @test body2.data == Dict(:meta => "data")

    atom1 = Atom(:Ar, [-2, -1, 3]u"bohr"; hello = :world)
    body3 = Atomistic.ElementMassBody(atom1)

    @test body3 isa Atomistic.ElementMassBody{Float64,Float64}
    @test body3.r ≈ @SVector [-2, -1, 3]
    @test body3.v ≈ @SVector [0, 0, 0]
    @test body3.m ≈ 39.9481austrip(1u"u")
    @test body3.symbol == :Ar
    @test body3.number == 18
    @test body3.data == Dict(:hello => :world)

    nbs_boundary_conditions = CubicPeriodicBoundaryConditions(1.5)
    atom2 = Atom(body3, nbs_boundary_conditions)

    @test position(atom2) ≈ (@SVector [1.0, 0.5, 0.0])u"bohr"
    @test velocity(atom2) ≈ velocity(atom1)
    @test atomic_mass(atom2) ≈ atomic_mass(atom1)
    @test AtomsBase.atomic_symbol(atom2) == AtomsBase.atomic_symbol(atom1)
    @test atomic_number(atom2) == atomic_number(atom1)
    @test atom2.data == atom1.data

    particles = [
        Atom(:Ar, [-2, -1, 3]u"bohr", [3, 5, 21]u"bohr * hartree / ħ_au"; a = :b),
        Atom(:Ar, [-7, 2, 13]u"bohr", [-3, 7, 1]u"bohr * hartree / ħ_au"; c = :d)
    ]
    box = (@SVector [(@SVector [1.5, 0.0, 0.0]), (@SVector [0.0, 1.5, 0.0]), (@SVector [0.0, 0.0, 1.5])])u"bohr"
    boundary_conditions = @SVector [Periodic(), Periodic(), Periodic()]
    system = FlexibleSystem(particles, box, boundary_conditions)
    bodies = Atomistic.get_bodies(system)
    flexible = FlexibleSystem(bodies, nbs_boundary_conditions)

    atoms_equal(a1::Atom, a2::Atom) = propertynames(a1) == propertynames(a2) && all(getproperty(a1, p) == getproperty(a2, p) for p ∈ propertynames(a1))
    bodies_equal(b1::Atomistic.ElementMassBody, b2::Atomistic.ElementMassBody) = all(getproperty(b1, p) == getproperty(b2, p) for p ∈ propertynames(b1))

    @test bodies_equal(bodies[1], Atomistic.ElementMassBody((@SVector [-2, -1, 3])u"bohr", (@SVector [3, 5, 21])u"bohr * hartree / ħ_au", elements[:Ar]; a = :b))
    @test bodies_equal(bodies[2], Atomistic.ElementMassBody((@SVector [-7, 2, 13])u"bohr", (@SVector [-3, 7, 1])u"bohr * hartree / ħ_au", elements[:Ar]; c = :d))

    @test atoms_equal(flexible[1], Atom(:Ar, [1.0, 0.5, 0.0]u"bohr", [3, 5, 21]u"bohr * hartree / ħ_au"; a = :b))
    @test atoms_equal(flexible[2], Atom(:Ar, [0.5, 0.5, 1.0]u"bohr", [-3, 7, 1]u"bohr * hartree / ħ_au"; c = :d))

    @test Atomistic.nbs_boundary_conditions(system) == nbs_boundary_conditions
    @test get_boundary_conditions(nbs_boundary_conditions) == boundary_conditions
    @test get_bounding_box(nbs_boundary_conditions) == box

    @test Atomistic.bound_position((@SVector [1, -1, 10]), nbs_boundary_conditions) == (@SVector [1, 0.5, 1])
end
