# Unit tests for implementations/molly/molly_simulator.jl

@testset "molly_simulator.jl" begin
    simulator1 = MollySimulator(400, 10, coupling=Molly.AndersenThermostat(94.4u"K", 0.1u"ps"))
    simulator2 = MollySimulator{Molly.StormerVerlet}(400, 10, t₀=1000, stride=2)
    simulator3 = MollySimulator{Molly.StormerVerlet}(400u"ns", 10)
    simulator4 = MollySimulator(400u"ns", 10, t₀=1000u"ns")

    @test simulator1 isa MollySimulator{Molly.VelocityVerlet,typeof(1u"ħ_au / hartree"),<:Molly.AndersenThermostat}
    @test simulator1.Δt == 400.0u"ħ_au / hartree"
    @test simulator1.steps == 10
    @test simulator1.t₀ == 0.0u"ħ_au / hartree"
    @test simulator1.coupling == Molly.AndersenThermostat(94.4u"K", 0.1u"ps")
    @test simulator1.stride == 1

    @test simulator2 isa MollySimulator{Molly.StormerVerlet,typeof(1u"ħ_au / hartree"),Molly.NoCoupling}
    @test simulator2.Δt == 400.0u"ħ_au / hartree"
    @test simulator2.steps == 10
    @test simulator2.t₀ == 1000.0u"ħ_au / hartree"
    @test simulator2.coupling == Molly.NoCoupling()
    @test simulator2.stride == 2

    @test simulator3 isa MollySimulator{Molly.StormerVerlet,typeof(1u"ns"),Molly.NoCoupling}
    @test simulator3.Δt == 400.0u"ns"
    @test simulator3.steps == 10
    @test simulator3.t₀ == 0.0u"ns"
    @test simulator3.coupling == Molly.NoCoupling()
    @test simulator3.stride == 1

    @test simulator4 isa MollySimulator{Molly.VelocityVerlet,typeof(1u"ns"),Molly.NoCoupling}
    @test simulator4.Δt == 400.0u"ns"
    @test simulator4.steps == 10
    @test simulator4.t₀ == 1000.0u"ns"
    @test simulator4.coupling == Molly.NoCoupling()
    @test simulator4.stride == 1

    particles = [
        AtomsBase.Atom(:Ar, (@SVector [7.0, 7.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 7.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 21.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [7.0, 21.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 7.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 7.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 21.0, 7.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au"),
        AtomsBase.Atom(:Ar, (@SVector [21.0, 21.0, 21.0])u"bohr", 6e-5(@SVector randn(3))u"bohr * hartree / ħ_au")
    ]
    box = [[28.0, 0.0, 0.0], [0.0, 28.0, 0.0], [0.0, 0.0, 28.0]]u"bohr"
    system = periodic_system(particles, box)

    potential = InteratomicPotentials.LennardJones(1.657e-21u"J", 0.34u"nm", 0.765u"nm", [:Ar])

    result = simulate(system, simulator1, potential)

    @test result isa MollyResult
end
