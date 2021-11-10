# Integrations with AtomsBase.jl

struct DynamicAtom{D,L <: Unitful.Length,V <: Unitful.Velocity} <: AbstractAtom
    position::SVector{D,L}
    velocity::SVector{D,V}
    element::ChemicalElement
end
function DynamicAtom(position, velocity, element)
    DynamicAtom{length(position),eltype(position),eltype(velocity)}(position, velocity, element)
end
function DynamicAtom(position, velocity, symbol::Union{Integer,AbstractString,Symbol,AbstractVector})
    DynamicAtom(position, velocity, ChemicalElement(symbol))
end

AtomsBase.position(atom::DynamicAtom) = atom.position
AtomsBase.velocity(atom::DynamicAtom) = atom.velocity
AtomsBase.element(atom::DynamicAtom) = atom.element
mass(atom::DynamicAtom) = atomic_mass(atom.element)

struct DynamicSystem{D,ET <: AbstractElement,AT <: AbstractParticle{ET},T <: Unitful.Length,TT <: Unitful.Time} <: AbstractSystem{D,ET,AT}
    box::SVector{D,SVector{D,T}}
    boundary_conditions::SVector{D,BoundaryCondition}
    particles::Vector{AT}
    time::TT
end

function DynamicSystem(box, boundary_conditions, particles, time)
    D = length(box)
    ET = typeof(element(first(particles)))
    AT = eltype(particles)
    T = eltype(first(box))
    TT = typeof(time)

    DynamicSystem{D,ET,AT,T,TT}(box, boundary_conditions, particles, time)
end

function DynamicSystem(particles::Vector{<:AbstractAtom}, box_size::Unitful.Length, time::Unitful.Time=0.0u"s")
    z = zero(typeof(box_size))
    box = SVector(SVector(box_size, z, z), SVector(z, box_size, z), SVector(z, z, box_size))
    boundary_conditions = SVector(Periodic(), Periodic(), Periodic())
    DynamicSystem(box, boundary_conditions, particles, time)
end

bounding_box(sys::DynamicSystem) = sys.box
boundary_conditions(sys::DynamicSystem) = sys.boundary_conditions

Base.size(sys::DynamicSystem) = size(sys.particles)
Base.length(sys::DynamicSystem) = length(sys.particles)
Base.getindex(sys::DynamicSystem, i::Int) = getindex(sys.particles, i)
