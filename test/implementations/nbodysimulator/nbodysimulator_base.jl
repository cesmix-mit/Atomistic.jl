# Unit tests for implementations/nbodysimulator/nbodysimulator_base.jl

@testset "implementations/nbodysimulator/nbodysimulator_base.jl" begin
    body1 = ElementMassBody((@SVector [1, 2, 3])u"nm", (@SVector [4, 5, 6])u"nm/s", :Ar; meta="data")

    @test body1 isa ElementMassBody{Float64,Float64}
    @test body1.r ≈ (@SVector [1, 2, 3])austrip(1u"nm")
    @test body1.v ≈ (@SVector [4, 5, 6])austrip(1u"nm/s")
    @test body1.m ≈ 39.9481austrip(1u"u")
    @test body1.symbol == :Ar
    @test body1.data == Dict(:meta => "data")

    atom1 = AtomsBase.Atom(:Ar, [-2, -1, 3]u"bohr"; hello=:world)
    body2 = ElementMassBody(atom1)

    @test body2 isa ElementMassBody{Float64,Float64}
    @test body2.r ≈ @SVector [-2, -1, 3]
    @test body2.v ≈ @SVector [0, 0, 0]
    @test body2.m ≈ 39.9481austrip(1u"u")
    @test body2.symbol == :Ar
    @test body2.data == Dict(:hello => :world)

    nbs_boundary_conditions = CubicPeriodicBoundaryConditions(1.5)
    atom2 = AtomsBase.Atom(body2, nbs_boundary_conditions)

    @test position(atom2) ≈ (@SVector [1.0, 0.5, 0.0])u"bohr"
    @test velocity(atom2) ≈ velocity(atom1)
    @test atomic_mass(atom2) ≈ atomic_mass(atom1)
    @test atomic_symbol(atom2) == atomic_symbol(atom1)
    @test atomic_number(atom2) == atomic_number(atom1)
    @test atom2.data == atom1.data

    atom3 = AtomsBase.Atom(body2, (@SVector [2, 3, 4]), (@SVector [5, 6, 7]), nbs_boundary_conditions)

    @test position(atom3) ≈ (@SVector [0.5, 0.0, 1.0])u"bohr"
    @test velocity(atom3) ≈ (@SVector [5, 6, 7])u"bohr * hartree / ħ_au"
    @test atomic_mass(atom3) ≈ atomic_mass(atom1)
    @test atomic_symbol(atom3) == atomic_symbol(atom1)
    @test atomic_number(atom3) == atomic_number(atom1)
    @test atom3.data == atom1.data

    atom4 = AtomsBase.Atom(body2, (@SVector [1.0, 0.5, 0.0])u"bohr", (@SVector [1.0, 0.5, 0.0])u"bohr/s")

    @test position(atom4) ≈ (@SVector [1.0, 0.5, 0.0])u"bohr"
    @test velocity(atom4) ≈ (@SVector [1.0, 0.5, 0.0])u"bohr/s"
    @test atomic_mass(atom4) ≈ atomic_mass(atom1)
    @test atomic_symbol(atom4) == atomic_symbol(atom1)
    @test atomic_number(atom4) == atomic_number(atom1)
    @test atom4.data == atom1.data

    particles = [
        AtomsBase.Atom(:Ar, [-2, -1, 3]u"bohr", [3, 5, 21]u"bohr * hartree / ħ_au"; a=:b),
        AtomsBase.Atom(:Ar, [-7, 2, 13]u"bohr", [-3, 7, 1]u"bohr * hartree / ħ_au"; c=:d)
    ]
    box = [[1.5, 0.0, 0.0], [0.0, 1.5, 0.0], [0.0, 0.0, 1.5]]u"bohr"
    system = periodic_system(particles, box)

    @test Atomistic.nbs_boundary_conditions(system) == nbs_boundary_conditions
    @test get_boundary_conditions(nbs_boundary_conditions) == [Periodic(), Periodic(), Periodic()]
    @test get_bounding_box(nbs_boundary_conditions) == box

    @test Atomistic.bound_position((@SVector [1, -1, 10]), nbs_boundary_conditions) == (@SVector [1, 0.5, 1])
end
