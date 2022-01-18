# Integrations with AtomsBase.jl

"""
    DynamicSystem{S<:AbstractSystem{D},T<:Unitful.Time} <: AbstractSystem{D}

Wraps any AbstractSystem with a time field to be able to cache expensive results.

**Type parameters**
- `D`: the dimension of the coordinate space
- `S`: the type of the wrapped system
- `T`: the type for time field
"""
struct DynamicSystem{D,S<:AbstractSystem{D},T<:Unitful.Time} <: AbstractSystem{D}
    system::S
    time::T
end

AtomsBase.bounding_box(s::DynamicSystem) = bounding_box(s.system)
AtomsBase.boundary_conditions(s::DynamicSystem) = boundary_conditions(s.system)
AtomsBase.species_type(s::DynamicSystem) = species_type(s.system)

Base.length(s::DynamicSystem) = length(s.system)
Base.getindex(s::DynamicSystem, i::Integer) = getindex(s.system, i)

AtomsBase.position(s::DynamicSystem) = position(s.system)
AtomsBase.velocity(s::DynamicSystem) = velocity(s.system)
AtomsBase.atomic_mass(s::DynamicSystem) = atomic_mass(s.system)
AtomsBase.atomic_symbol(s::DynamicSystem) = AtomsBase.atomic_symbol(s.system)
AtomsBase.atomic_number(s::DynamicSystem) = atomic_number(s.system)

AtomsBase.position(s::DynamicSystem, i) = position(s.system, i)
AtomsBase.velocity(s::DynamicSystem, i) = velocity(s.system, i)
AtomsBase.atomic_mass(s::DynamicSystem, i) = atomic_mass(s.system, i)
AtomsBase.atomic_symbol(s::DynamicSystem, i) = AtomsBase.atomic_symbol(s.system, i)
AtomsBase.atomic_number(s::DynamicSystem, i) = atomic_number(s.system, i)
