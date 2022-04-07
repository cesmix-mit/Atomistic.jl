# Unit tests for implementations/molly/molly_result.jl

@testset "molly_result.jl" begin
    particles = [
        AtomsBase.Atom(:Ar, (@SVector [7.0, 7.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"; meta=:data),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 7.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"; hello="world"),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 21.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 21.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 7.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 7.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 21.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 21.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au")
    ]
    box = (@SVector [(@SVector [28.0, 0.0, 0.0]), (@SVector [0.0, 28.0, 0.0]), (@SVector [0.0, 0.0, 28.0])])u"bohr"
    boundary_conditions = @SVector [Periodic(), Periodic(), Periodic()]
    system = FlexibleSystem(particles, box, boundary_conditions)

    simulator = MollySimulator(400, 10, t₀=1000, coupling=Molly.AndersenThermostat(94.4u"K", 0.1u"ps"))
    simulator2 = MollySimulator(400, 9, t₀=1000, stride=2)
    simulator3 = MollySimulator(400, 10, t₀=1000, stride=2)

    potential = InteratomicPotentials.LennardJones(austrip(1.657e-21u"J"), austrip(0.34u"nm"), austrip(0.765u"nm"), [:Ar])

    result = simulate(system, simulator, potential)
    result2 = simulate(system, simulator2, potential)
    result3 = simulate(system, simulator3, potential)

    @test get_time_range(result) == (1000:400:4600)u"ħ_au / hartree"
    @test get_time_range(result2) == (1400:800:3800)u"ħ_au / hartree"
    @test get_time_range(result3) == (1400:800:4600)u"ħ_au / hartree"

    @test get_num_bodies(result) == 8
    @test get_bounding_box(result) == box
    @test get_boundary_conditions(result) == boundary_conditions

    @test reference_temperature(result) == 94.4u"K"
    @test ismissing(reference_temperature(result2))

    @test get_positions(result) isa AbstractVector{<:StaticVector{3,<:Unitful.Length}}
    @test get_velocities(result) isa AbstractVector{<:StaticVector{3,<:Unitful.Velocity}}
    @test get_particles(result) isa AbstractVector{<:AtomsBase.Atom}
    @test get_system(result) isa System{3}

    @test all(all(0.0u"bohr" ≤ c < 28.0u"bohr" for c ∈ p) for p ∈ get_positions(result))
    @test get_particles(result)[1].meta == :data && get_particles(result)[2].hello == "world"

    @test all(get_positions(result, t) isa AbstractVector{<:StaticVector{3,<:Unitful.Length}} for t ∈ 1:10)
    @test all(get_velocities(result, t) isa AbstractVector{<:StaticVector{3,<:Unitful.Velocity}} for t ∈ 1:10)
    @test all(get_particles(result, t) isa AbstractVector{<:AtomsBase.Atom} for t ∈ 1:10)

    @test all(all(all(0.0u"bohr" ≤ c < 28.0u"bohr" for c ∈ p) for p ∈ get_positions(result, t)) for t ∈ 1:10)
    @test all(get_particles(result, t)[1].meta == :data && get_particles(result, t)[2].hello == "world" for t ∈ 1:10)

    @test Atomistic.temperature(result) isa Unitful.Temperature
    @test Atomistic.kinetic_energy(result) isa Unitful.Energy
    @test Atomistic.potential_energy(result) isa Unitful.Energy

    @test all(isapprox(Atomistic.total_energy(result3, t), Atomistic.total_energy(result3), rtol=0.1) for t ∈ 1:5)
end
