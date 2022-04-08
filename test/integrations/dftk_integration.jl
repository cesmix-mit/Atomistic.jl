# Unit tests for integrations/dftk_integration.jl

@testset "integrations/dftk_integration.jl" begin
    potential = DFTKPotential(5u"hartree", [1, 1, 1]; functionals=[:lda_x, :lda_c_pw], scf_kwargs=Dict(:damping => 0.7, :tol => 1e-4))

    particles = [AtomsBase.Atom(:Ar, [i & 1 == 0 ? 7 : 21, i & 2 == 0 ? 7 : 21, i & 4 == 0 ? 7 : 21]u"bohr") for i ∈ 0:7]
    box = [[28.0, 0.0, 0.0], [0.0, 28.0, 0.0], [0.0, 0.0, 28.0]]u"bohr"
    system = attach_psp(periodic_system(particles, box); functional="lda")

    eandf = energy_and_force(system, potential)
    @test eandf.e isa Float64
    @test eandf.f isa AbstractVector{<:SVector{3,<:Float64}}

    @test haskey(potential.scf_kwargs, :ψ)
    @test haskey(potential.scf_kwargs, :ρ)
end
