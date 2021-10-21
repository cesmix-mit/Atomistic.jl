# Integrations with DFTK.jl
# This integration should ultimately live within the DFTK package itself

@kwdef struct DFTKPotential <: ArbitraryPotential
    psp::ElementPsp
    lattice
    Ecut::Real
    kgrid::AbstractVector{Integer}
    n_bands::Union{Integer,Nothing} = nothing
    tol::Union{AbstractFloat,Nothing} = nothing
    damping::Union{AbstractFloat,Nothing} = nothing
    mixing::Union{Mixing,Nothing} = nothing
    previous_scfres::RefValue{Any} = Ref{Any}()
end
function DFTKPotential(;
                       psp::ElementPsp,
                       lattice,
                       Ecut::Quantity,
                       kgrid::AbstractVector{Integer},
                       n_bands::Union{Integer,Nothing}=nothing,
                       tol::Union{AbstractFloat,Nothing}=nothing,
                       damping::Union{AbstractFloat,Nothing}=nothing,
                       mixing::Union{Mixing,Nothing}=nothing,
                       previous_scfres::RefValue{Any}=Ref{Any}())
    DFTKPotential(psp, lattice, austrip(Ecut), kgrid, n_bands, tol, damping, mixing, previous_scfres)
end

function InteratomicPotentials.force(state::MassBodies, potential::DFTKPotential)
    compute_forces_cart(calculate_scf(state, potential))[1]
end

function potential_energy(state::MassBodies, potential::DFTKPotential)
    calculate_scf(state, potential).energies.total
end

function calculate_scf(state::MassBodies, potential::DFTKPotential)
    state = DFTKAtoms(state, potential.psp, potential.lattice)
    model = model_LDA(state.lattice, state.atoms)
    basis = PlaneWaveBasis(model; Ecut=potential.Ecut, kgrid=potential.kgrid)

    extra_args = isassigned(potential.previous_scfres) ? (ψ = potential.previous_scfres[].ψ, ρ = potential.previous_scfres[].ρ) : (; )
    scfres = self_consistent_field(basis; extra_args..., (f => getfield(potential, f) for f ∈ (:n_bands, :tol, :damping, :mixing) if getfield(potential, f) !== nothing)...)
    potential.previous_scfres[] = scfres
end

function analyze_convergence(state::MassBodies, potential::DFTKPotential, cutoffs::Vector{<:Quantity})
    energies = Vector{Float64}()
    for Ecut ∈ cutoffs
        parameters = DFTKPotential(psp=potential.psp,
                                   lattice=potential.lattice,
                                   Ecut=Ecut,
                                   kgrid=potential.kgrid,
                                   n_bands=potential.n_bands,
                                   tol=potential.tol,
                                   damping=potential.damping,
                                   mixing=potential.mixing)
        scfres = calculate_scf(state, parameters)
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
