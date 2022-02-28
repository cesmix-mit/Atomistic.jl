# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# TODO: all the auconverts can be removed once InteratomicPotentials.jl properly supports unitful
# Convert arbitrary AbstractSystem to a System
function Molly.System(system::AbstractSystem{D}; kwargs...) where {D}
    @assert hcat(bounding_box(system)...) == bounding_box(system)[1][1] * I(D)
    @assert all(periodicity(system))
    atoms_data = (species_type(system) <: AtomsBase.Atom) ? [AugmentedAtomData(AtomsBase.atomic_symbol(a); a.data...) for a ∈ system] :
                 (species_type(system) <: Molly.Atom) ? copy(system.atoms_data) :
                 AugmentedAtomData.(AtomsBase.atomic_symbol(system))
    velocities = ismissing(AtomsBase.velocity(system)) ? [@SVector zeros(VELOCITY_TYPE, 3) for _ ∈ 1:length(system)] :
                 [auconvert.(v) for v ∈ AtomsBase.velocity(system)]
    System(;
        atoms = [Molly.Atom(index = i, mass = auconvert(m)) for (i, m) ∈ enumerate(atomic_mass(system))],
        atoms_data = atoms_data,
        coords = [auconvert.(p) for p ∈ position(system)],
        velocities = velocities,
        box_size = SVector{D}(auconvert(bounding_box(system)[i][i]) for i ∈ 1:D),
        force_units = FORCE_UNIT,
        energy_units = ENERGY_UNIT,
        kwargs...
    )
end

# Alternative atom data type to store arbitrary fields
struct AugmentedAtomData
    element::Symbol
    data::Dict{Symbol,Any}
end
AugmentedAtomData(element::Symbol; data...) = AugmentedAtomData(element, Dict{Symbol,Any}(data...))
Base.hasproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) || haskey(data.data, x)
Base.getproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) ? getfield(data, x) : getindex(data.data, x)
function Base.propertynames(data::AugmentedAtomData, private::Bool = false)
    if private
        (fieldnames(AugmentedAtomData)..., keys(data.data)...)
    else
        (filter(!isequal(:data), fieldnames(AugmentedAtomData))..., keys(data.data)...)
    end
end

# Convert Molly Atom to AtomsBase Atom
AtomsBase.Atom(d::AugmentedAtomData, p::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}) = AtomsBase.Atom(d.element, p, v; d.data...)

AtomsBase.atomic_symbol(s::System, i) = Symbol(s.atoms_data[i].element)
