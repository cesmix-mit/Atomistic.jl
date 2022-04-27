# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# Convert arbitrary AbstractSystem to a System
function Molly.System(system::AbstractSystem{D}; kwargs...) where {D}
    @assert hcat(bounding_box(system)...) == bounding_box(system)[1][1] * I(D)
    @assert all(periodicity(system))
    atoms_data = (species_type(system) <: AtomsBase.Atom) ? [AugmentedAtomData(atomic_symbol(a); a.data...) for a ∈ system] :
                 (species_type(system) <: Molly.Atom) ? copy(system.atoms_data) :
                 AugmentedAtomData.(atomic_symbol(system))
    velocities = ismissing(velocity(system)) ? [@SVector zeros(VELOCITY_TYPE, 3) for _ ∈ 1:length(system)] :
                 velocity(system)
    System(;
        atoms=[Molly.Atom(index=i, mass=m) for (i, m) ∈ enumerate(atomic_mass(system))],
        atoms_data,
        coords=position(system),
        velocities,
        box_size=SVector{D}(bounding_box(system)[i][i] for i ∈ 1:D),
        force_units=FORCE_UNIT,
        energy_units=ENERGY_UNIT,
        kwargs...
    )
end

# Alternative atom data type to store arbitrary fields
struct AugmentedAtomData
    element::Symbol
    data::Dict{Symbol,Any}
end
AugmentedAtomData(element::Symbol; data...) = AugmentedAtomData(element, Dict(data...))
Base.hasproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) || haskey(data.data, x)
Base.getproperty(data::AugmentedAtomData, x::Symbol) = hasfield(AugmentedAtomData, x) ? getfield(data, x) : getindex(data.data, x)
function Base.propertynames(data::AugmentedAtomData, private::Bool=false)
    if private
        (fieldnames(AugmentedAtomData)..., keys(data.data)...)
    else
        (filter(!isequal(:data), fieldnames(AugmentedAtomData))..., keys(data.data)...)
    end
end

Base.hasproperty(view::AtomView{<:System}, x::Symbol) = hasfield(AtomView, x) || haskey(getfield(view, :system).atoms_data[getfield(view, :index)].data, x)
Base.getproperty(view::AtomView{<:System}, x::Symbol) = hasfield(AtomView, x) ? getfield(view, x) : getindex(getfield(view, :system).atoms_data[getfield(view, :index)].data, x)
Base.propertynames(view::AtomView{<:System}, private::Bool=false) = (fieldnames(AtomView)..., keys(getfield(view, :system).atoms_data[getfield(view, :index)].data)...)

# Convert Molly Atom to AtomsBase Atom
AtomsBase.Atom(d::AugmentedAtomData, p::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}) = AtomsBase.Atom(d.element, p, v; d.data...)
