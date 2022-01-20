# Molecular Dynamics Result Abstract Interface

"""
    MolecularDynamicsResult
Abstract type to be extended by all concrete structs representing the result of a molecular dynamics simulation.
"""
abstract type MolecularDynamicsResult end

"""
    get_system(result::MolecularDynamicsResult, t::Integer = 0)::AbstractSystem

Extract the underlying system at a particular timestep from the simulation result.

An implementer of this API should implement a method of this function for their custom result type.
"""
function get_system(result::MolecularDynamicsResult, t::Integer = 0)
    throw(UnimplementedError(:get_system, result))
end
"""
    get_time_range(result::MolecularDynamicsResult)::AbstractVector{Unitful.Time}

Extract the unit-anotated time range from the simulation result.

An implementer of this API should implement a method of this function for their custom result type.
"""
function get_time_range(result::MolecularDynamicsResult)
    throw(UnimplementedError(:get_time_range, result))
end

"""
    temperature(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Temperature

Extract the unit-anotated temperature of the simulation at a particular timestep from the simulation result.
The default sentinel value of `t = 0` indicates the _end_ of the simulation.

An implementer of this API should implement a method of this function for their custom result type.
"""
function temperature(result::MolecularDynamicsResult, t::Integer = 0)
    throw(UnimplementedError(:temperature, result))
end
"""
    reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}

Extract the unit-anotated reference temperature of the simulation from the simulation result.
If there is no thermostat with a reference temperature in this simulation, return missing.

An implementer of this API should implement a method of this function for their custom result type if it supports thermostats.
If not implmented, the default implemenation just returns missing.
"""
function reference_temperature(result::MolecularDynamicsResult)
    missing
end

"""
    kinetic_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy

Extract the unit-anotated kinetic energy of the simulation at a particular timestep from the simulation result.
The default sentinel value of `t = 0` indicates the _end_ of the simulation.

An implementer of this API should implement a method of this function for their custom result type.
"""
function kinetic_energy(result::MolecularDynamicsResult, t::Integer = 0)
    throw(UnimplementedError(:kinetic_energy, result))
end
"""
    potential_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy

Extract the unit-anotated potential energy of the simulation at a particular timestep from the simulation result.
The default sentinel value of `t = 0` indicates the _end_ of the simulation.

An implementer of this API should implement a method of this function for their custom result type.
"""
function potential_energy(result::MolecularDynamicsResult, t::Integer = 0)
    throw(UnimplementedError(:potential_energy, result))
end
"""
    total_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy

Extract the unit-anotated total energy of the simulation at a particular timestep from the simulation result.
The default sentinel value of `t = 0` indicates the _end_ of the simulation.

The default implementation simply sums the kinetic and potential energy at the timestep.
An implementer of this API could implement a method of this function for their custom result type if it supports a more efficient way to calculate this quantity.
"""
function total_energy(result::MolecularDynamicsResult, t::Integer = 0)
    kinetic_energy(result, t) + potential_energy(result, t)
end

"""
    rdf(result::MolecularDynamicsResult, sample_fraction::Float64 = 1.0)::Tuple{AbstractVector{Real},AbstractVector{Real}}

Calculate the radial distribution function from the simulation result.

To include only a trailing portion of the timesteps for reduced noise and faster computation, set `sample_fraction` to be less than 1.0; `sample_fraction` must be in the range (0.0, 1.0].
The result is a tuple of vectors which represent the interparticle radial distances (in bohr) and the density of each distance respectively.
"""
function rdf(result::MolecularDynamicsResult, sample_fraction::Float64 = 1.0)
    throw(UnimplementedError(:rdf, result))
end