# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

function Molly.System(system::AbstractSystem{D}; kwargs...) where {D}
    @assert hcat(bounding_box(system)...) == bounding_box(system)[1][1] * I(D)
    @assert all(periodicity(system))
    atoms_data = (species_type(system) <: AtomsBase.Atom) ? [AugmentedAtomData(AtomsBase.atomic_symbol(a), a.data) for a ∈ system] :
                 (species_type(system) <: Molly.Atom) ? copy(system.atoms_data) :
                 [AugmentedAtomData(AtomsBase.atomic_symbol(a), Dict{Symbol,Any}()) for a ∈ system]
    System(;
        atoms = [Molly.Atom(index = i, mass = atomic_mass(a)) for (i, a) ∈ enumerate(system)],
        atoms_data = atoms_data,
        coords = position(system),
        velocities = AtomsBase.velocity(system),
        box_size = SVector{D}([bounding_box(system)[i][i] for i ∈ 1:D]),
        force_units = FORCE_UNIT,
        energy_units = ENERGY_UNIT,
        kwargs...
    )
end

struct AugmentedAtomData
    element::Symbol
    data::Dict{Symbol,Any}
end

AtomsBase.Atom(d::AtomData, p::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}) = AtomsBase.Atom(d.element, p, v; data...)

AtomsBase.atomic_symbol(s::System, i) = Symbol(s.atoms_data[i].element)
