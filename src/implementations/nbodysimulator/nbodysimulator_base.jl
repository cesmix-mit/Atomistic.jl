# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# Internal type used to represent a MassBody and retain the element metadata
struct ElementMassBody{cType<:Real,mType<:Real} <: Body
    r::SVector{3,cType} # in LENGTH_UNIT
    v::SVector{3,cType} # in VELOCITY_UNIT
    m::mType            # in MASS_UNIT
    symbol::Symbol
    data::Dict{Symbol,Any}
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, symbol::Union{Integer,AbstractString,Symbol}; data...)
    r, v = promote(austrip.(r), austrip.(v))
    e = elements[symbol]
    ElementMassBody(r, v, austrip(e.atomic_mass), Symbol(e.symbol), Dict{Symbol,Any}(data...))
end

# Convert AtomsBase Atom to NBodySimulator body
function ElementMassBody(atom::AtomsBase.Atom)
    r, v = promote(austrip.(position(atom)), austrip.(AtomsBase.velocity(atom)))
    ElementMassBody(r, v, austrip(atomic_mass(atom)), AtomsBase.atomic_symbol(atom), atom.data)
end
# Convert NBodySimulator body to AtomsBase Atom
AtomsBase.Atom(b::ElementMassBody, boundary_conditions::BoundaryConditions) = AtomsBase.Atom(b, b.r, b.v, boundary_conditions)
AtomsBase.Atom(b::ElementMassBody, r::SVector{3,<:Real}, v::SVector{3,<:Real}, boundary_conditions::BoundaryConditions) = AtomsBase.Atom(b, bound_position(r, boundary_conditions) * LENGTH_UNIT, v * VELOCITY_UNIT)
AtomsBase.Atom(b::ElementMassBody, r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}) = AtomsBase.Atom(b.symbol, r, v; b.data...)

# Convert AtomsBase boundary conditions to NBodySimulator boundary conditions
function nbs_boundary_conditions(system::AbstractSystem{3})
    # TODO: support more boundary conditions
    @assert hcat(bounding_box(system)...) == bounding_box(system)[1][1] * I(3)
    @assert all(periodicity(system))
    CubicPeriodicBoundaryConditions(austrip(bounding_box(system)[1][1]))
end
# Convert NBodySimulator boundary conditions to AtomsBase boundary conditions
# TODO: support more boundary conditions
get_boundary_conditions(::CubicPeriodicBoundaryConditions) = @SVector [Periodic(), Periodic(), Periodic()]
# Convert NBodySimulator boundary conditions to AtomsBase bounding box
# TODO: support more boundary conditions
function get_bounding_box(boundary_conditions::CubicPeriodicBoundaryConditions)
    SVector{3}(SVector{3}.(eachrow(boundary_conditions.L * LENGTH_UNIT * I(3))))
end

# Bound a position according to the boundary conditions
# TODO: support more boundary conditions
bound_position(r::SVector{3,<:Real}, boundary_conditions::CubicPeriodicBoundaryConditions) = mod.(r, boundary_conditions.L)
