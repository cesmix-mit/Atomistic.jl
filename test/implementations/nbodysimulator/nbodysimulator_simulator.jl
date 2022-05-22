# Unit tests for implementations/nbodysimulator/nbodysimulator_simulator.jl

@testset "nbodysimulator_simulator.jl" begin
    simulator1 = NBSimulator(400, 10, thermostat=NBodySimulator.AndersenThermostat(100, 2e-4))
    simulator2 = NBSimulator(400, 10, t₀=1000)
    simulator3 = NBSimulator(400u"ns", 10, simulator=DPRKN6())
    simulator4 = NBSimulator(400u"ns", 10, t₀=1000u"ns")

    @test simulator1 isa NBSimulator{NBodySimulator.VelocityVerlet,typeof(1u"ħ_au / hartree"),<:NBodySimulator.AndersenThermostat}
    @test simulator1.Δt == 400.0u"ħ_au / hartree"
    @test simulator1.steps == 10
    @test simulator1.t₀ == 0.0u"ħ_au / hartree"
    @test simulator1.thermostat == NBodySimulator.AndersenThermostat(100, 2e-4)
    @test simulator1.simulator == NBodySimulator.VelocityVerlet()
    @test Atomistic.time_range(simulator1) == (0.0, 4000.0)

    @test simulator2 isa NBSimulator{NBodySimulator.VelocityVerlet,typeof(1u"ħ_au / hartree"),NBodySimulator.NullThermostat}
    @test simulator2.Δt == 400.0u"ħ_au / hartree"
    @test simulator2.steps == 10
    @test simulator2.t₀ == 1000.0u"ħ_au / hartree"
    @test simulator2.thermostat == NBodySimulator.NullThermostat()
    @test simulator2.simulator == NBodySimulator.VelocityVerlet()
    @test Atomistic.time_range(simulator2) == (1000.0, 5000.0)

    @test simulator3 isa NBSimulator{DPRKN6,typeof(1u"ns"),NBodySimulator.NullThermostat}
    @test simulator3.Δt == 400.0u"ns"
    @test simulator3.steps == 10
    @test simulator3.t₀ == 0.0u"ns"
    @test simulator3.thermostat == NBodySimulator.NullThermostat()
    @test simulator3.simulator == DPRKN6()
    @test all(Atomistic.time_range(simulator3) .≈ (0.0, 10austrip(400.0u"ns")))

    @test simulator4 isa NBSimulator{NBodySimulator.VelocityVerlet,typeof(1u"ns"),NBodySimulator.NullThermostat}
    @test simulator4.Δt == 400.0u"ns"
    @test simulator4.steps == 10
    @test simulator4.t₀ == 1000.0u"ns"
    @test simulator4.thermostat == NBodySimulator.NullThermostat()
    @test simulator4.simulator == NBodySimulator.VelocityVerlet()
    @test all(Atomistic.time_range(simulator4) .≈ (austrip(1000.0u"ns"), austrip(1000.0u"ns") + 10austrip(400.0u"ns")))

    particles = [
        AtomsBase.Atom(:Ar, (@SVector [7, 7, 7])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7, 7, 21])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7, 21, 7])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7, 21, 21])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21, 7, 7])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21, 7, 21])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21, 21, 7])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21, 21, 21])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au")
    ]
    box = [[28.0, 0.0, 0.0], [0.0, 28.0, 0.0], [0.0, 0.0, 28.0]]u"bohr"
    system = periodic_system(particles, box)

    potential = InteratomicPotentials.LennardJones(1.657e-21u"J", 0.34u"nm", 0.765u"nm", [:Ar])

    result = simulate(system, simulator1, potential)

    @test result isa NBSResult
    @test length(result.energy_cache) == length(result.result.solution.t)
end
