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
        atoms = [Molly.Atom(index = i, mass = auconvert(atomic_mass(a))) for (i, a) ∈ enumerate(system)],
        atoms_data = atoms_data,
        coords = [auconvert.(p) for p ∈ position(system)],
        velocities = [auconvert.(v) for v ∈ AtomsBase.velocity(system)],
        box_size = SVector{D}(auconvert(bounding_box(system)[i][i]) for i ∈ 1:D),
        force_units = FORCE_UNIT,
        energy_units = ENERGY_UNIT,
        kwargs...
    )
end

struct AugmentedAtomData
    element::Symbol
    data::Dict{Symbol,Any}
end
Base.hasproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) || haskey(data.data, x)
Base.getproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) ? getfield(data, x) : getindex(data.data, x)
function Base.propertynames(data::AugmentedAtomData, private::Bool = false)
    if private
        (fieldnames(AugmentedAtomData)..., keys(data.data)...)
    else
        (filter(!isequal(:data), fieldnames(AugmentedAtomData))..., keys(data.data)...)
    end
end

AtomsBase.Atom(d::AugmentedAtomData, p::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}) = AtomsBase.Atom(d.element, p, v; d.data...)

AtomsBase.atomic_symbol(s::System, i) = Symbol(s.atoms_data[i].element)
