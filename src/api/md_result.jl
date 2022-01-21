# Molecular Dynamics Result Abstract Interface

"""
    MolecularDynamicsResult
Abstract type to be extended by all concrete structs representing the result of a molecular dynamics simulation.
"""
abstract type MolecularDynamicsResult end

# -----------------------------------------------------------------------------
# Simulation Configuration
# -----------------------------------------------------------------------------

"""
    get_time_range(result::MolecularDynamicsResult)::AbstractVector{Unitful.Time}

Extract the unit-anotated time range of the simulation from the simulation result.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_time_range(result::MolecularDynamicsResult) = throw(UnimplementedError(:get_time_range, result))
Base.length(result::MolecularDynamicsResult) = length(get_time_range(result::MolecularDynamicsResult))

"""
    get_bounding_box(result::MolecularDynamicsResult)::SVector{3, SVector{3, Unitful.Length}}

Extract the unit-anotated bounding box of the simulation from the simulation result in an `AtomBase`-compatible format.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_bounding_box(result::MolecularDynamicsResult) = throw(UnimplementedError(:get_bounding_box, result))

"""
    get_boundary_conditions(result::MolecularDynamicsResult)::SVector{3, BoundaryCondition}

Extract the boundary conditions of the simulation from the simulation result in an `AtomBase`-compatible format.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_boundary_conditions(result::MolecularDynamicsResult) = throw(UnimplementedError(:get_boundary_conditions, result))

"""
    reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}

Extract the unit-anotated reference temperature of the simulation from the simulation result.
If there is no thermostat with a reference temperature in this simulation, return missing.

An implementer of this API should implement a method of this function for their custom result type if it supports thermostats.
If not implmented, the default implemenation just returns missing.
"""
reference_temperature(result::MolecularDynamicsResult) = missing

# -----------------------------------------------------------------------------
# Time-Dependent System Data
# -----------------------------------------------------------------------------

"""
    get_time(result::MolecularDynamicsResult, t::Integer)::Unitful.Time

Extract the unit-anotated time of the simulation at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

The default implementation extracts the time from the result of `get_time_range`.
"""
get_time(result::MolecularDynamicsResult, t::Integer) = get_time_range(result)[t]
get_time(result::MolecularDynamicsResult) = get_time(result, length(result))

"""
    get_positions(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Length}}

Extract the unit-anotated position of each particle in the system at a particular timestep from the simulation result in an `AtomBase`-compatible format.
The returned positions should be normalized such that they lie within the bounding box of the system, even if the system is periodic.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_positions(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:get_positions, result))
get_positions(result::MolecularDynamicsResult) = get_positions(result, length(result))

"""
    get_velocities(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Velocity}}

Extract the unit-anotated velocity of each particle in the system at a particular timestep from the simulation result in an `AtomBase`-compatible format.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_velocities(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:get_velocities, result))
get_velocities(result::MolecularDynamicsResult) = get_velocities(result, length(result))

"""
    get_particles(result::MolecularDynamicsResult, t::Integer)::AbstractVector{Atom}

Extract the position of each particle in the system at a particular timestep from the simulation result in an `AtomBase`-compatible format.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
get_particles(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:get_particles, result))
get_particles(result::MolecularDynamicsResult) = get_particles(result, length(result))

"""
    get_system(result::MolecularDynamicsResult, t::Integer)::AbstractSystem{3}

Extract the underlying system at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

The default implementation combines the results of `get_particles`, `get_bounding_box`, `get_boundary_conditions`, and `get_time` at the timestep to create a `FlexibleSystem` wrapped in a `DyanmicSystem`.
An implementer of this API could implement a method of this function for their custom result type if it supports a more efficient way to extract the system.
"""
function get_system(result::MolecularDynamicsResult, t::Integer)
    system = FlexibleSystem(get_particles(result, t), get_bounding_box(result), get_boundary_conditions(result))
    DynamicSystem(system, get_time(result, t))
end
get_system(result::MolecularDynamicsResult) = get_system(result, length(result))

# -----------------------------------------------------------------------------
# Simulation Analysis
# -----------------------------------------------------------------------------

"""
    temperature(result::MolecularDynamicsResult, t::Integer)::Unitful.Temperature

Extract the unit-anotated temperature of the simulation at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
temperature(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:temperature, result))
temperature(result::MolecularDynamicsResult) = temperature(result, length(result))

"""
    kinetic_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy

Extract the unit-anotated kinetic energy of the simulation at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
kinetic_energy(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:kinetic_energy, result))
kinetic_energy(result::MolecularDynamicsResult) = kinetic_energy(result, length(result))

"""
    potential_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy

Extract the unit-anotated potential energy of the simulation at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

An implementer of this API should implement a method of this function for their custom result type.
"""
potential_energy(result::MolecularDynamicsResult, t::Integer) = throw(UnimplementedError(:potential_energy, result))
potential_energy(result::MolecularDynamicsResult) = potential_energy(result, length(result))

"""
    total_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy

Extract the unit-anotated total energy of the simulation at a particular timestep from the simulation result.
Automatically defaults to the end of the simulation when `t` is not passed.

The default implementation simply sums the results of `kinetic_energy` and `potential_energy` at the timestep.
An implementer of this API could implement a method of this function for their custom result type if it supports a more efficient way to calculate this quantity.
"""
total_energy(result::MolecularDynamicsResult, t::Integer) = kinetic_energy(result, t) + potential_energy(result, t)
total_energy(result::MolecularDynamicsResult) = total_energy(result, length(result))

"""
    rdf(result::MolecularDynamicsResult, sample_fraction::Float64 = 1.0)::Tuple{AbstractVector{Real},AbstractVector{Real}}

Calculate the radial distribution function from the simulation result.

To include only a trailing portion of the timesteps for reduced noise and faster computation, set `sample_fraction` to be less than 1.0; `sample_fraction` must be in the range (0.0, 1.0].
The result is a tuple of vectors which represent the interparticle radial distances (in bohr) and the density of each distance respectively.
"""
rdf(result::MolecularDynamicsResult, sample_fraction::Float64 = 1.0) = throw(UnimplementedError(:rdf, result))
