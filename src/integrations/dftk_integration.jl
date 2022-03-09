# Integrations with DFTK.jl

# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# ! all functions in this section are copied directly from experimental DFTK code
# ! this code will be deleted once it is live in DFTK

function parse_system(system::AbstractSystem{3})
    if !all(periodicity(system))
        error("DFTK only supports calculations with periodic boundary conditions.")
    end

    # Parse abstract system and return data required to construct model
    lattice = austrip.(hcat(bounding_box(system)...))
    T = eltype(lattice)

    # Cache for instantiated pseudopotentials
    # (such that the respective objects are indistinguishable)
    cached_pseudos = Dict{String,Any}()
    atoms = map(system) do atom
        if hasproperty(atom, :potential)
            potential = atom.potential
        elseif hasproperty(atom, :pseudopotential)
            pspkey = atom.pseudopotential
            potential = get!(cached_pseudos, pspkey) do
                ElementPsp(AtomsBase.atomic_symbol(atom); psp = load_psp(pspkey))
            end
        else
            potential = ElementCoulomb(AtomsBase.atomic_symbol(atom))
        end

        potential => SVector{3,T}(lattice \ T.(austrip.(position(atom))))
    end

    oldatoms = oldatoms_from_new(atoms)

    (; lattice, atoms = oldatoms)
end
function oldatoms_from_new(atomic_potentials)
    potentials = first.(atomic_potentials)
    potential_groups = [findall(Ref(pot) .== potentials) for pot in Set(potentials)]
    [first(atomic_potentials[first(group)]) => last.(atomic_potentials[group]) for group in potential_groups]
end
function DFTK.model_LDA(system::AbstractSystem; kwargs...)
    parsed = parse_system(system)
    model_LDA(parsed.lattice, parsed.atoms; kwargs...)
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

# ! This integration should ultimately move to the DFTK package itself

@kwdef struct DFTKPotential <: ArbitraryPotential
    Ecut::Real
    kgrid::AbstractVector{<:Integer}
    n_bands::Union{Integer,Nothing} = nothing
    tol::Union{AbstractFloat,Nothing} = nothing
    damping::Union{AbstractFloat,Nothing} = nothing
    mixing::Union{Mixing,Nothing} = nothing
    previous_scfres::Ref{Any} = Ref{Any}()
    potential_energy_cache::Dict{Float64,Float64} = Dict{Float64,Float64}()
end
function DFTKPotential(Ecut::Real, kgrid::AbstractVector{<:Integer}; kwargs...)
    DFTKPotential(; Ecut = Ecut, kgrid = kgrid, kwargs...)
end
function DFTKPotential(Ecut::Unitful.Energy, kgrid::AbstractVector{<:Integer}; kwargs...)
    DFTKPotential(; Ecut = austrip(Ecut), kgrid = kgrid, kwargs...)
end

function InteratomicPotentials.energy_and_force(system::AbstractSystem, potential::DFTKPotential)
    model = model_LDA(system)
    basis = PlaneWaveBasis(model; Ecut = potential.Ecut, kgrid = potential.kgrid)

    args = (f => getfield(potential, f) for f ∈ (:n_bands, :tol, :damping, :mixing) if !isnothing(getfield(potential, f)))
    extra_args = isassigned(potential.previous_scfres) ? (ψ = potential.previous_scfres[].ψ, ρ = potential.previous_scfres[].ρ) : (;)
    scfres = self_consistent_field(basis; args..., extra_args...)
    potential.previous_scfres[] = scfres
    # TODO: support multiple species
    (; e = scfres.energies.total, f = compute_forces_cart(scfres)[1])
end
