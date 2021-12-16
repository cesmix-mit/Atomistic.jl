# Integrations with AtomsBase.jl

abstract type AbstractAtom end

"""
    DynamicAtom{D,L<:Unitful.Length,V<:Unitful.Velocity} <: AbstractAtom

An atom representation based on the StaticAtom but with a velocity.

**Type parameters**
- `D`: the dimension of the coordinate space
- `L`: the type for the position components
- `V`: the type for the velocity components

Fields should not be accessed directly. Use the provided accessors instead.
"""
struct DynamicAtom{D,L<:Unitful.Length,V<:Unitful.Velocity} <: AbstractAtom
    position::SVector{D,L}
    velocity::SVector{D,V}
    element::Element
end
function DynamicAtom(position, velocity, element)
    DynamicAtom{length(position),eltype(position),eltype(velocity)}(position, velocity, element)
end
function DynamicAtom(position, velocity, symbol::Union{Integer,AbstractString,Symbol,AbstractVector})
    DynamicAtom(position, velocity, elements(symbol))
end

AtomsBase.position(atom::DynamicAtom) = atom.position
AtomsBase.velocity(atom::DynamicAtom) = atom.velocity
AtomsBase.species(atom::DynamicAtom) = atom.element

AtomsBase.atomic_symbol(a::DynamicAtom) = a.element.symbol
AtomsBase.atomic_mass(a::DynamicAtom) = a.element.atomic_mass
AtomsBase.atomic_number(a::DynamicAtom) = a.element.number
AtomsBase.atomic_property(a::DynamicAtom, property::Symbol) = getproperty(a.element, property)

"""
    DynamicSystem{D,AT<:AbstractAtom,L<:Unitful.Length,TT<:Unitful.Time} <: AbstractAtomicSystem{D}

A representation of a system of dynamic atoms which is similar to the FlexibleSystem but with a time field.

**Type parameters**
- `D`: the dimension of the coordinate space
- `A`: the type for the atoms that make up the system
- `L`: the type for the bounding box components
- `T`: the type for time field

Fields should not be accessed directly. Use the provided accessors instead.
"""
struct DynamicSystem{D,A<:AbstractAtom,L<:Unitful.Length,T<:Unitful.Time} <: AbstractAtomicSystem{D}
    box::SVector{D,SVector{D,L}}
    boundary_conditions::SVector{D,BoundaryCondition}
    particles::Vector{A}
    time::T
end
function DynamicSystem(box, boundary_conditions, particles, time)
    D = length(box)
    A = eltype(particles)
    L = eltype(eltype(box))
    T = typeof(time)

    DynamicSystem{D,A,L,T}(box, boundary_conditions, particles, time)
end
function DynamicSystem(particles::Vector{<:AbstractAtom}, box_size::Unitful.Length, time::Unitful.Time = 0.0u"s")
    z = zero(typeof(box_size))
    box = SVector(SVector(box_size, z, z), SVector(z, box_size, z), SVector(z, z, box_size))
    boundary_conditions = SVector(Periodic(), Periodic(), Periodic())
    DynamicSystem(box, boundary_conditions, particles, time)
end

AtomsBase.bounding_box(sys::DynamicSystem) = sys.box
AtomsBase.boundary_conditions(sys::DynamicSystem) = sys.boundary_conditions

Base.size(sys::DynamicSystem) = size(sys.particles)
Base.length(sys::DynamicSystem) = length(sys.particles)
Base.getindex(sys::DynamicSystem, i::Int) = getindex(sys.particles, i)
