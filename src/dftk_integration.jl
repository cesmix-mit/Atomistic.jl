setup_threading(n_blas=4)

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
    previous_scfres::Base.RefValue{Any} = Ref{Any}()
end
DFTKForceGenerationParameters(parameters::DFTKForceGenerationParameters, Ecut::Quantity) = DFTKForceGenerationParameters(
    box_size=parameters.box_size,
    psp=parameters.psp,
    lattice=parameters.lattice,
    Ecut=Ecut,
    kgrid=parameters.kgrid,
    n_bands=parameters.n_bands,
    tol=parameters.tol,
    α=parameters.α,
    mixing=parameters.mixing,
)

function dftk_atoms(element::DFTK.Element, bodies::AbstractVector{MassBody}, box_size::Quantity)
    [element => [austrip.(b.r * u"bohr" / box_size) for b ∈ bodies]]
end

function calculate_scf(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters)
    atoms = dftk_atoms(parameters.psp, bodies, parameters.box_size)

    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    extra_args = isassigned(parameters.previous_scfres) ? (ψ=parameters.previous_scfres[].ψ, ρ=parameters.previous_scfres[].ρ) : (; )
    scfres =  @time self_consistent_field(basis; extra_args..., (f=>getfield(parameters, f) for f ∈ (:n_bands, :tol, :α, :mixing) if getfield(parameters, f) !== missing)...)
    parameters.previous_scfres[] = scfres
    return scfres
end

function generate_forces(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters)
    scfres = calculate_scf(bodies, parameters)
    return compute_forces_cart(scfres)[1]
end

function analyze_convergence(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters, cutoffs::AbstractVector{<:Quantity})
    energies = Vector{Float64}()
    for Ecut ∈ cutoffs
        params = DFTKForceGenerationParameters(parameters, Ecut)
        scfres = calculate_scf(bodies, params)
        push!(energies, scfres.energies.total)
    end
    
    plot(
		title="DFTK Analysis",
		xlab="Ecut",
		ylab="Total Energy",
        legend=false,
        cutoffs,
        energies * u"hartree"
	)
end
