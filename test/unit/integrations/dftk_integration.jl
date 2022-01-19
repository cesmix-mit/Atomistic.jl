# Unit tests for integrations/dftk_integration.jl

@testset "integrations/dftk_integration.jl" begin
    pspkey = list_psp(:Ar, functional = "lda")[1].identifier
    particles = [Atom(:Ar, [i & 1 == 0 ? 7 : 21, i & 2 == 0 ? 7 : 21, i & 4 == 0 ? 7 : 21]u"bohr"; pseudopotential = pspkey) for i ∈ 0:7]
    box = [[28.0, 0, 0], [0, 28.0, 0], [0, 0, 28.0]]u"bohr"
    boundary_conditions = [Periodic(), Periodic(), Periodic()]
    system1 = FlexibleSystem(particles, box, boundary_conditions)

    potential1 = DFTKPotential(5, [1, 1, 1]; damping = 0.7, tol = 1e-4)
    potential2 = DFTKPotential(136.05693123044495u"eV", [1, 1, 1]; damping = 0.7, tol = 1e-4)

    @test potential1.Ecut ≈ potential2.Ecut

    @test Atomistic.calculate_scf(system1, potential1) == potential1.previous_scfres[]
    @test InteratomicPotentials.potential_energy(system1, potential1) isa Float64
    @test InteratomicPotentials.force(system1, potential1) isa AbstractVector{<:StaticVector{3,<:Float64}}

    system2 = DynamicSystem(system1, 0u"s")
    system3 = DynamicSystem(system1, 1u"ħ_au / hartree")

    @test InteratomicPotentials.potential_energy(system2, potential1) == potential1.potential_energy_cache[0]
    @test InteratomicPotentials.force(system3, potential1) isa AbstractVector{<:StaticVector{3,<:Float64}}
    @test keys(potential1.potential_energy_cache) == Set([0.0, 1.0])
    pe2 = potential1.potential_energy_cache[0]
    pe3 = potential1.potential_energy_cache[1]

    # workaround way to validate that cached results are used
    modified_particles = [Atom(:Ar, [i & 1 == 0 ? 5 : 15, i & 2 == 0 ? 5 : 15, i & 4 == 0 ? 5 : 15]u"bohr"; pseudopotential = pspkey) for i ∈ 0:7]
    system4 = FlexibleSystem(modified_particles, box, boundary_conditions)
    system5 = DynamicSystem(system4, 0u"s")
    system6 = DynamicSystem(system4, 1u"ħ_au / hartree")

    @test InteratomicPotentials.potential_energy(system5, potential1) == pe2
    @test InteratomicPotentials.potential_energy(system6, potential1) == pe3
    @test keys(potential1.potential_energy_cache) == Set([0.0, 1.0])
end
