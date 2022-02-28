# Unit tests for integrations/dftk_integration.jl

@testset "integrations/dftk_integration.jl" begin
    potential1 = DFTKPotential(5, [1, 1, 1]; damping = 0.7, tol = 1e-4)
    potential2 = DFTKPotential(136.05693123044495u"eV", [1, 1, 1]; damping = 0.7, tol = 1e-4)

    @test potential1.Ecut ≈ potential2.Ecut

    pspkey = list_psp(:Ar, functional = "lda")[1].identifier
    particles = [AtomsBase.Atom(:Ar, [i & 1 == 0 ? 7 : 21, i & 2 == 0 ? 7 : 21, i & 4 == 0 ? 7 : 21]u"bohr"; pseudopotential = pspkey) for i ∈ 0:7]
    box = [[28.0, 0, 0], [0, 28.0, 0], [0, 0, 28.0]]u"bohr"
    boundary_conditions = [Periodic(), Periodic(), Periodic()]
    system1 = FlexibleSystem(particles, box, boundary_conditions)
    eandf1 = energy_and_force(system1, potential1)

    @test eandf1.e isa Float64
    @test eandf1.f isa AbstractVector{<:SVector{3,<:Float64}}

    @test potential1.previous_scfres[].energies.total == eandf1.e
    @test compute_forces_cart(potential1.previous_scfres[])[1] == eandf1.f

    system2 = System(system1)
    eandf2 = energy_and_force(system2, potential2)

    @test eandf2.e isa Float64
    @test eandf2.f isa AbstractVector{<:SVector{3,<:Float64}}

    @test potential2.previous_scfres[].energies.total == eandf2.e
    @test compute_forces_cart(potential2.previous_scfres[])[1] == eandf2.f
end
