# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# Internal type used to represent a MassBody and retain the element metadata
struct ElementMassBody{cType<:Real,mType<:Real} <: Body
    r::SVector{3,cType} # in LENGTH_UNIT
    v::SVector{3,cType} # in VELOCITY_UNIT
    m::mType            # in MASS_UNIT
    s::Symbol
    n::Int
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, e::Element)
    ElementMassBody{Float64,Float64}(austrip.(r), austrip.(v), austrip(e.atomic_mass), Symbol(e.symbol), e.number)
end

# Convert AtomsBase Atom to NBodySimulator body
function ElementMassBody(atom::Atom)
    ElementMassBody{Float64,Float64}(austrip.(position(atom)), austrip.(velocity(atom)), austrip(atomic_mass(atom)), AtomsBase.atomic_symbol(atom), atomic_number(atom))
end
# Convert NBodySimulator body to AtomsBase Atom
# TODO: support more boundary conditions
function AtomsBase.Atom(b::ElementMassBody, boundary_conditions::CubicPeriodicBoundaryConditions)
    Atom(b.s, mod.(b.r, boundary_conditions.L) .* LENGTH_UNIT, b.v .* VELOCITY_UNIT)
end

# Convert AtomsBase boundary conditions to NBodySimulator boundary conditions
function nbs_boundary_conditions(system::AbstractSystem)
    # TODO: support more boundary conditions
    box = bounding_box(system)
    @assert box[1][1] == box[2][2] == box[3][3]
    @assert boundary_conditions(system) == [Periodic(), Periodic(), Periodic()]
    CubicPeriodicBoundaryConditions(austrip(box[1][1]))
end
# Convert NBodySimulator boundary conditions to AtomsBase boundary conditions
# TODO: support more boundary conditions
ab_boundary_conditions(::CubicPeriodicBoundaryConditions) = [Periodic(), Periodic(), Periodic()]

# Convert NBodySimulator boundary conditions to AtomsBase bounding box
# TODO: support more boundary conditions
function ab_bounding_box(boundary_conditions::CubicPeriodicBoundaryConditions)
    box_size = boundary_conditions.L * LENGTH_UNIT
    z = zero(typeof(box_size))
    [[box_size, z, z], [z, box_size, z], [z, z, box_size]]
end

# Convert AtomsBase AbstractSystem to Vector of NBodySimulator bodies
function bodies(system::AbstractSystem)
    ElementMassBody.(system)
end
# Convert Vector of NBodySimulator bodies to AtomsBase FlexibleSystem
function AtomsBase.FlexibleSystem(bodies::AbstractVector{<:ElementMassBody}, boundary_conditions::BoundaryConditions)
    particles = Fix2(Atom, boundary_conditions).(bodies)
    FlexibleSystem(particles, ab_bounding_box(boundary_conditions), ab_boundary_conditions(boundary_conditions))
end
# Convert Vector of NBodySimulator bodies to AtomsBase FastSystem
function AtomsBase.FastSystem(bodies::AbstractVector{<:ElementMassBody}, boundary_conditions::BoundaryConditions)
    particles = Fix2(Atom, boundary_conditions).(bodies)
    FastSystem(particles, ab_bounding_box(boundary_conditions), ab_boundary_conditions(boundary_conditions))
end

# Extract the atomic symbol of each body in the system
get_atomic_symbols(system::PotentialNBodySystem{ElementMassBody}) = [b.s for b ∈ system.bodies]
# Extract the atomic number of each body in the system
get_atomic_numbers(system::PotentialNBodySystem{ElementMassBody}) = [b.n for b ∈ system.bodies]

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
