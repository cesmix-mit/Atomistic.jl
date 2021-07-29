# Integrations with DFTK.jl

Base.@kwdef struct DFTKParameters <: NuclearPotentialParameters
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
    psp=parameters.psp,
    lattice=parameters.lattice,
    Ecut=Ecut,
    kgrid=parameters.kgrid,
    n_bands=parameters.n_bands,
    tol=parameters.tol,
    α=parameters.α,
    mixing=parameters.mixing,
)

function forces(state::AtomicConfiguration, parameters::DFTKParameters)
    compute_forces_cart(calculate_scf(state, parameters))[1]
end

function potential_energy(state::AtomicConfiguration, parameters::DFTKParameters)
    calculate_scf(state, parameters).energies.total
end

function calculate_scf(state::AtomicConfiguration, parameters::DFTKParameters)
    setup_threading(n_blas=4)
    
    state = DFTKAtoms(state, parameters.psp, parameters.lattice)
    model = model_LDA(state.lattice, state.atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    extra_args = isassigned(parameters.previous_scfres) ? (ψ=parameters.previous_scfres[].ψ, ρ=parameters.previous_scfres[].ρ) : (; )
    scfres = self_consistent_field(basis; extra_args..., (f=>getfield(parameters, f) for f ∈ (:n_bands, :tol, :α, :mixing) if getfield(parameters, f) !== nothing)...)
    parameters.previous_scfres[] = scfres
end

function analyze_convergence(state::AtomicConfiguration, parameters::DFTKParameters, cutoffs::AbstractVector{<:Quantity})
    energies = Vector{Float64}()
    for Ecut ∈ cutoffs
        params = DFTKParameters(parameters, Ecut)
        scfres = calculate_scf(state, params)
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
