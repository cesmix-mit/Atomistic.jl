# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# Internal type used to represent a MassBody and retain the element metadata
struct ElementMassBody{cType<:Real,mType<:Real} <: Body
    r::SVector{3,cType} # in LENGTH_UNIT
    v::SVector{3,cType} # in VELOCITY_UNIT
    m::mType            # in MASS_UNIT
    symbol::Symbol
    number::Int
    data::Dict{Symbol,Any}
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, e::Element; data...)
    r, v = promote(austrip.(r), austrip.(v))
    ElementMassBody(r, v, austrip(e.atomic_mass), Symbol(e.symbol), e.number, Dict{Symbol,Any}(data...))
end
function ElementMassBody(body::ElementMassBody, r::SVector{3,<:Real}, v::SVector{3,<:Real})
    r, v = promote(r, v)
    ElementMassBody(r, v, body.m, body.symbol, body.number, body.data)
end

# Convert AtomsBase Atom to NBodySimulator body
function ElementMassBody(atom::Atom)
    r, v = promote(austrip.(position(atom)), austrip.(velocity(atom)))
    ElementMassBody(r, v, austrip(atomic_mass(atom)), AtomsBase.atomic_symbol(atom), atomic_number(atom), atom.data)
end
# Convert NBodySimulator body to AtomsBase Atom
function AtomsBase.Atom(b::ElementMassBody, boundary_conditions::CubicPeriodicBoundaryConditions)
    Atom(b.symbol, bound_position(b.r, boundary_conditions) .* LENGTH_UNIT, b.v .* VELOCITY_UNIT; b.data...)
end

# Convert AtomsBase AbstractSystem to Vector of NBodySimulator bodies
get_bodies(system::AbstractSystem{3}) = ElementMassBody.(system)
# Convert Vector of NBodySimulator bodies to AtomsBase FlexibleSystem
function AtomsBase.FlexibleSystem(bodies::AbstractVector{<:ElementMassBody}, boundary_conditions::BoundaryConditions)
    particles = Fix2(Atom, boundary_conditions).(bodies)
    FlexibleSystem(particles, get_bounding_box(boundary_conditions), get_boundary_conditions(boundary_conditions))
end

# Convert AtomsBase boundary conditions to NBodySimulator boundary conditions
function nbs_boundary_conditions(system::AbstractSystem{3})
    # TODO: support more boundary conditions
    @assert hcat(bounding_box(system)...) == Matrix(I * bounding_box(system)[1][1], 3, 3)
    @assert all(periodicity(system))
    CubicPeriodicBoundaryConditions(austrip(bounding_box(system)[1][1]))
end
# Convert NBodySimulator boundary conditions to AtomsBase boundary conditions
# TODO: support more boundary conditions
get_boundary_conditions(::CubicPeriodicBoundaryConditions) = @SVector [Periodic(), Periodic(), Periodic()]
# Convert NBodySimulator boundary conditions to AtomsBase bounding box
# TODO: support more boundary conditions
function get_bounding_box(boundary_conditions::CubicPeriodicBoundaryConditions)
    SVector{3}(SVector{3}.(eachrow(Matrix(I * (boundary_conditions.L * LENGTH_UNIT), 3, 3))))
end

# Bound a position according to the boundary conditions
# TODO: support more boundary conditions
bound_position(r::SVector{3,<:Real}, boundary_conditions::CubicPeriodicBoundaryConditions) = mod.(r, boundary_conditions.L)


# -----------------------------------------------------------------------------
# Convenience methods for generating starting configurations with element data
# -----------------------------------------------------------------------------

function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, symbol::Union{Integer,AbstractString,Symbol}, L::Unitful.Length, reference_temp::Unitful.Temperature; rng = MersenneTwister(n))
    e = elements[symbol]
    average_velocity = √(u"k" * reference_temp / e.atomic_mass)
    generate_bodies_in_cell_nodes(n, symbol, average_velocity, L, rng = rng)
end
function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, symbol::Union{Integer,AbstractString,Symbol}, average_velocity::Unitful.Velocity, L::Unitful.Length; rng = MersenneTwister(n))
    velocities = average_velocity * randn(rng, Float64, (3, n))
    e = elements[symbol]
    bodies = ElementMassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x ∈ dL/2:dL:L, y ∈ dL/2:dL:L, z ∈ dL/2:dL:L
        if count > n
            break
        end
        push!(bodies, ElementMassBody(SVector(x, y, z), SVector{3}(velocities[:, count]), e))
        count += 1
    end
    return bodies
end
