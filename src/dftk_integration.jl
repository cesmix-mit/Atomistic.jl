setup_threading()

Base.@kwdef struct DFTKForceGenerationParameters <: ForceGenerationParameters
    box_size::Quantity
    psp::ElementPsp
    lattice::AbstractArray{Quantity, 2}
    Ecut::Quantity
    kgrid::AbstractVector{Integer}
    n_bands::Union{Integer, Missing} = missing
    tol::Union{AbstractFloat, Missing} = missing
    α::Union{AbstractFloat, Missing} = missing
    mixing = missing # There is no abstract type for mixing :(
end
DFTKForceGenerationParameters(parameters::DFTKForceGenerationParameters, Ecut::Quantity) = DFTKForceGenerationParameters(parameters.box_size, parameters.psp, parameters.lattice, Ecut, parameters.kgrid, parameters.n_bands, parameters.tol, parameters.α, parameters.mixing);

function calculate_scf(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters)
    atoms = [parameters.psp => [auconvert.(u"bohr", b.r) / parameters.box_size for b ∈ bodies]]
    
    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    return @time self_consistent_field(basis; (f=>getfield(parameters, f) for f in (:n_bands, :tol, :α, :mixing) if getfield(parameters, f) !== missing)...)
end

function generate_forces(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters)
    scfres = calculate_scf(bodies, parameters)
    return ParticleForcePotentialParameters(compute_forces_cart(scfres)[1])
end

function analyze_convergence(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters, cutoffs::Vector{<:Quantity})
    options = Dict(Ecut => DFTKForceGenerationParameters(parameters, Ecut) for Ecut in cutoffs)
    fields = Dict(Ecut => calculate_scf(bodies, options[Ecut]) for Ecut in cutoffs)
    energies = Dict(Ecut => fields[Ecut].energies.total for Ecut in cutoffs)
    plot(
		title="DFTK Analysis",
		xlab="Ecut",
		ylab="Total Energy",
        legend=false
	)
	plot!(
        cutoffs,
        c -> auconvert(u"hartree", energies[c])
	)
end
