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

function AtomsBase.SimpleSystem(particles::Vector{<:AbstractAtom}, box_size::Unitful.Length)
    z = zero(typeof(box_size))
    box = SVector(SVector(box_size, z, z), SVector(z, box_size, z), SVector(z, z, box_size))
    boundary_conditions = SVector(Periodic(), Periodic(), Periodic())
    SimpleSystem(box, boundary_conditions, particles)
end
