# Unit tests for implementations/nbodysimulator/nbodysimulator_simulator.jl

@testset "nbodysimulator_simulator.jl" begin
    system = FlexibleSystem([
            Atomistic.ElementMassBody((@SVector [0, 0, 1])u"bohr", (@SVector [0, 0, -1])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, 1, 0])u"bohr", (@SVector [0, -1, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [1, 0, 0])u"bohr", (@SVector [-1, 0, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, 0, -1])u"bohr", (@SVector [0, 0, 1])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [0, -1, 0])u"bohr", (@SVector [0, 1, 0])u"bohr * hartree / ħ_au", elements[:Ar]),
            Atomistic.ElementMassBody((@SVector [-1, 0, 0])u"bohr", (@SVector [1, 0, 0])u"bohr * hartree / ħ_au", elements[:Ar])
        ], CubicPeriodicBoundaryConditions(60.0))

    simulator1 = NBSimulator(400, 10, thermostat = AndersenThermostat(100, 0.5))
    simulator2 = NBSimulator(400, 10, t₀ = 1000)
    simulator3 = NBSimulator(400u"ns", 10, simulator = DPRKN6())
    simulator4 = NBSimulator(400u"ns", 10, t₀ = 1000u"ns")

    potential1 = LennardJonesParameters(1.657e-21u"J", 0.34u"nm", 0.765u"nm")
    potential2 = LennardJones(austrip(1.657e-21u"J"), austrip(0.34u"nm"), austrip(0.765u"nm"))

    @test simulator1.Δt == 400.0
    @test simulator1.steps == 10
    @test simulator1.t₀ == 0.0
    @test simulator1.thermostat == AndersenThermostat(100, 0.5)
    @test simulator1.simulator == VelocityVerlet()

    @test simulator2.Δt == 400.0
    @test simulator2.steps == 10
    @test simulator2.t₀ == 1000.0
    @test simulator2.thermostat == NBodySimulator.NullThermostat()
    @test simulator2.simulator == VelocityVerlet()

    @test simulator3.Δt == austrip(400.0u"ns")
    @test simulator3.steps == 10
    @test simulator3.t₀ == 0.0
    @test simulator3.thermostat == NBodySimulator.NullThermostat()
    @test simulator3.simulator == DPRKN6()

    @test simulator4.Δt == austrip(400.0u"ns")
    @test simulator4.steps == 10
    @test simulator4.t₀ == austrip(1000.0u"ns")
    @test simulator4.thermostat == NBodySimulator.NullThermostat()
    @test simulator4.simulator == VelocityVerlet()

    @test simulate(system, simulator1, potential1) isa NBSResult
    @test simulate(system, simulator2, potential2) isa NBSResult
end
