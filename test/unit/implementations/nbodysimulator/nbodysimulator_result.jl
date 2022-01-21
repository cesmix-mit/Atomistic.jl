# Unit tests for implementations/nbodysimulator/nbodysimulator_result.jl

@testset "nbodysimulator_result.jl" begin
    system = FlexibleSystem([
            Atomistic.ElementMassBody((@SVector [0, 0, 1])u"bohr", (@SVector [0, 0, -1])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, 1, 0])u"bohr", (@SVector [0, -1, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [1, 0, 0])u"bohr", (@SVector [-1, 0, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, 0, -1])u"bohr", (@SVector [0, 0, 1])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, -1, 0])u"bohr", (@SVector [0, 1, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [-1, 0, 0])u"bohr", (@SVector [1, 0, 0])u"bohr * hartree / ħ_au", elements[:Ar])
        ], CubicPeriodicBoundaryConditions(60.0))

    simulator1 = NBSimulator(400, 10, t₀ = 1000, thermostat = AndersenThermostat(0.0002, 0.0002))
    simulator2 = NBSimulator(400, 10, t₀ = 1000)

    potential1 = LennardJonesParameters(1.657e-21u"J", 0.34u"nm", 0.765u"nm")
    potential2 = LennardJones(austrip(1.657e-21u"J"), austrip(0.34u"nm"), austrip(0.765u"nm"))

    result = simulate(system, simulator1, potential1)
    result1 = simulate(system, simulator2, potential1)
    result2 = simulate(system, simulator2, potential2)

    @test get_time_range(result1) == (1000:400:5000)u"ħ_au / hartree"
    @test get_bounding_box(result1) == (@SVector [(@SVector [60.0, 0.0, 0.0]), (@SVector [0.0, 60.0, 0.0]), (@SVector [0.0, 0.0, 60.0])])u"bohr"
    @test get_boundary_conditions(result1) == @SVector [Periodic(), Periodic(), Periodic()]

    @test reference_temperature(result) == 0.0002u"hartree / k_au"
    @test ismissing(reference_temperature(result1))

    @test get_positions(result) isa AbstractVector{<:StaticVector{3,<:Unitful.Length}}
    @test all(all(all(0u"bohr" ≤ c < 60u"bohr" for c ∈ p) for p ∈ get_positions(result, t)) for t ∈ 1:11)
    @test get_velocities(result) isa AbstractVector{<:StaticVector{3,<:Unitful.Velocity}}
    @test get_particles(result) isa AbstractVector{<:Atom}

    @test Atomistic.temperature(result) isa Unitful.Temperature
    @test Atomistic.kinetic_energy(result) isa Unitful.Energy
    @test Atomistic.potential_energy(result) isa Unitful.Energy

    # @test all(Atomistic.total_energy(result1, t) ≈ Atomistic.total_energy(result1) for t ∈ 1:10)
    # @test all(Atomistic.total_energy(result2, t) ≈ Atomistic.total_energy(result2) for t ∈ 1:10)
end
