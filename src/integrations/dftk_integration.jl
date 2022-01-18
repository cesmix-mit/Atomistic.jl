# Integrations with DFTK.jl

# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# ! all functions in this section are copied directly from experimental DFTK code
# ! this code will be deleted once it is live in DFTK

function parse_system(system::AbstractSystem{D}) where {D}
    if !all(periodicity(system))
        error("DFTK only supports calculations with periodic boundary conditions.")
    end

    # Parse abstract system and return data required to construct model
    mtx = austrip.(hcat(bounding_box(system)...))
    T = eltype(mtx)
    lattice = zeros(T, 3, 3)
    lattice[1:D, 1:D] .= mtx

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

        coordinate = zeros(T, 3)
        coordinate[1:D] = inv(lattice[1:D, 1:D]) * T.(austrip.(position(atom)))
        potential => Vec3{T}(coordinate)
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
    previous_scfres::RefValue{Any} = Ref{Any}()
    potential_energy_cache::Dict{Float64,Float64} = Dict{Float64,Float64}()
end
function DFTKPotential(
    Ecut::Unitful.Energy,
    kgrid::AbstractVector{<:Integer},
    n_bands::Union{Integer,Nothing} = nothing,
    tol::Union{AbstractFloat,Nothing} = nothing,
    damping::Union{AbstractFloat,Nothing} = nothing,
    mixing::Union{Mixing,Nothing} = nothing,
    previous_scfres::RefValue{Any} = Ref{Any}(),
    potential_energy_cache::Dict{Float64,Float64} = Dict{Float64,Float64}()
)
    DFTKPotential(austrip(Ecut), kgrid, n_bands, tol, damping, mixing, previous_scfres, potential_energy_cache)
end

function calculate_scf(system::AbstractSystem, potential::DFTKPotential)
    model = model_LDA(system)
    basis = PlaneWaveBasis(model; Ecut = potential.Ecut, kgrid = potential.kgrid)

    args = (f => getfield(potential, f) for f ∈ (:n_bands, :tol, :damping, :mixing) if !isnothing(getfield(potential, f)))
    extra_args = isassigned(potential.previous_scfres) ? (ψ = potential.previous_scfres[].ψ, ρ = potential.previous_scfres[].ρ) : (;)
    scfres = self_consistent_field(basis; args..., extra_args...)
    potential.previous_scfres[] = scfres
end

function InteratomicPotentials.potential_energy(system::AbstractSystem, potential::DFTKPotential)
    calculate_scf(system, potential).energies.total
end

function InteratomicPotentials.potential_energy(system::DynamicSystem, potential::DFTKPotential)
    get!(potential.potential_energy_cache, austrip(system.time)) do
        calculate_scf(system, potential).energies.total
    end
end

function InteratomicPotentials.force(system::AbstractSystem, potential::DFTKPotential)
    # TODO: support multiple species
    compute_forces_cart(calculate_scf(system, potential))[1]
end

function InteratomicPotentials.force(system::DynamicSystem, potential::DFTKPotential)
    scf = calculate_scf(system, potential)
    potential.potential_energy_cache[austrip(system.time)] = scf.energies.total
    # TODO: support multiple species
    compute_forces_cart(scf)[1]
end

# -----------------------------------------------------------------------------
# Experimental functions
# -----------------------------------------------------------------------------

# ! This function is experimental and will eventually be removed

function analyze_convergence(system::AbstractSystem, potential::DFTKPotential, cutoffs::AbstractVector{<:Unitful.Energy})
    energies = zeros(Float64, length(cutoffs))
    for i ∈ 1:length(cutoffs)
        parameters = DFTKPotential(
            Ecut = cutoffs[i],
            kgrid = potential.kgrid,
            n_bands = potential.n_bands,
            tol = potential.tol,
            damping = potential.damping,
            mixing = potential.mixing
        )
        @info "Ecut: $(cutoffs[i])"
        energies[i] = InteratomicPotentials.potential_energy(system, parameters)
    end

    plot(
        title = "DFTK Analysis",
        xlab = "Ecut",
        ylab = "Total Energy",
        legend = false,
        cutoffs,
        energies * ENERGY_UNIT
    )
end
