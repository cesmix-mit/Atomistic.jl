# Integrations with DFTK.jl

Base.@kwdef struct DFTKParameters <: NuclearPotentialParameters
    box_size::Quantity
    psp::ElementPsp
    lattice::AbstractArray{Quantity, 2}
    Ecut::Quantity
    kgrid::AbstractVector{Integer}
    n_bands::Union{Integer, Nothing} = nothing
    tol::Union{AbstractFloat, Nothing} = nothing
    α::Union{AbstractFloat, Nothing} = nothing
    mixing = nothing # There is no abstract type for mixing :(
    previous_scfres::Base.RefValue{Any} = Ref{Any}()
end
DFTKParameters(parameters::DFTKParameters, Ecut::Quantity) = DFTKParameters(
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

function generate_forces(bodies::AbstractVector{<:MassBody}, parameters::DFTKParameters)
    scfres = calculate_scf(bodies, parameters)
    return compute_forces_cart(scfres)[1]
end

function calculate_scf(bodies::AbstractVector{<:MassBody}, parameters::DFTKParameters)
    setup_threading(n_blas=4)
    
    atoms = dftk_atoms(parameters.psp, bodies, parameters.box_size)

    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    extra_args = isassigned(parameters.previous_scfres) ? (ψ=parameters.previous_scfres[].ψ, ρ=parameters.previous_scfres[].ρ) : (; )
    scfres =  @time self_consistent_field(basis; extra_args..., (f=>getfield(parameters, f) for f ∈ (:n_bands, :tol, :α, :mixing) if getfield(parameters, f) !== nothing)...)
    parameters.previous_scfres[] = scfres
    return scfres
end

function dftk_atoms(element::DFTK.Element, bodies::AbstractVector{<:MassBody}, box_size::Quantity)
    [element => [austrip.(b.r * u"bohr" / box_size) for b ∈ bodies]]
end

function analyze_convergence(bodies::AbstractVector{<:MassBody}, parameters::DFTKParameters, cutoffs::AbstractVector{<:Quantity})
    energies = Vector{Float64}()
    for Ecut ∈ cutoffs
        params = DFTKParameters(parameters, Ecut)
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
