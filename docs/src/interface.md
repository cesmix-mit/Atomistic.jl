# Implementing the Atomistic API

There are two main parts to the Atomistic API. An implementer must provide a concrete simulator type subtyping `MolecularDynamicsSimulator` and a concrete result type subtyping `MolecularDynamicsResult`. A summary of the required functionality for each component follows.

# Molecular Dynamics Simulators

The only function that must be implemented for subtypes of `MolecularDynamicsSimulator` is a method of `simulate`.

```julia
simulate(system::AbstractSystem, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult
```

The first parameter is an instance of an implementation of the `AbstractSystem` interface defined by [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl). It is important that the implementer rely only on the general API rather than any specific implementation details because the user could provide any system implementation; if there are performance-sensitive decisions based on the system type, this should be achieved through multiple dispatch. The species, position, and velocity of each body as well as the boundary shape and conditions are stored in this data structure. An implementer should take special care to handle units appropriately when converting from this representation to the internal representation of the underlying simulator package.

The second parameter is the custom simulator. This could specify information such as the duration of the simulation, numerical details such as discretization parameters and integration methods, or domain specific options such as thermostats and barostats.

The third parameter is an instance of an implementation of the `ArbitraryPotential` interface defined by [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl). The main function from this interface that is useful in the dynamics context is `InteratomicPotentials.force(s::AbstractSystem, p::ArbitraryPotential)::AbstractVector{StaticVector{3, Real}}`. However, it is recommended that implementers consider using `InteratomicPotentials.energy_and_force(s::AbstractSystem, p::ArbitraryPotential)::NamedTuple{(:e, :f), Tuple{Real, Vector{SVector{3, Real}}}}` and cache the energy results on the `MolecularDynamicsResult` struct because the additional energy calculation is asymptotically free for most interatomic potential implementations. See the discussion of the `potential_energy` function below.

The `simulate` method should return the corresponding implementation of the `MolecularDynamicsResult` interface as described below.

# [Molecular Dynamics Results](@id MolecularDynamicsResult_Specification)

Fifteen total functions must be implemented for subtypes of `MolecularDynamicsResult`, four of which have default implementations. The functions allow users to access simulation data from three categories: simulation configuration, time-dependent system data, and simulation analysis.

## Simulation Configuration Functions

There are five functions in the `MolecularDynamicsResult` API which return time-independent simulation configuration information:

```julia
get_time_range(result::MolecularDynamicsResult)::AbstractVector{Unitful.Time}
```

This function should return a unit-anotated vector containing the time value for each step of the simulation which can be iterated for plotting, animation, or other analysis.

```julia
get_num_bodies(result::MolecularDynamicsResult)::Integer
```

This function should return the numebr of bodies in the simulation.

```julia
get_bounding_box(result::MolecularDynamicsResult)::SVector{3, SVector{3, Unitful.Length}}
```

This function should return the [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl)-style bounding box of the system used in the simulation.

```julia
get_boundary_conditions(result::MolecularDynamicsResult)::SVector{3, BoundaryCondition}
```

This function should return the [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl)-style boundary conditions of the system used in the simulation.

```julia
reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}
```

This function should return a unit-anotated temperature that describes the temperature maintained by the thermostat used in the simulation. If there is no reference temperature for the particular simulation, the function should return missing, which is the default implementation.

## Time-Dependent System Data

There are five functions in the `MolecularDynamicsResult` API which return system data at a particular timestep, two of which have default implementations.

```julia
get_time(result::MolecularDynamicsResult, t::Integer)::Unitful.Time
```

This convenience function converts a timestep to the unit-annotated time of the timestep using the provided `get_time_range` implementation. The timestep defaults to the end of the simulation when `t` is not passed. Most implementers will not need to write a custom implementation.

```julia
get_positions(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Length}}
```

This function should return the unit-annotated position of each particle in the system at a particular timestep from the simulation result. The returned positions should be normalized such that they lie within the bounding box of the system, even if the system is periodic. The timestep defaults to the end of the simulation when `t` is not passed.

```julia
get_velocities(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Velocity}}
```

This function should return the unit-annotated velocities of each particle in the system at a particular timestep from the simulation result. The timestep defaults to the end of the simulation when `t` is not passed.

```julia
get_particles(result::MolecularDynamicsResult, t::Integer)::AbstractVector{Atom}
```

This function should return the particles in the system at a particular timestep from the simulation result. The timestep defaults to the end of the simulation when `t` is not passed. It is important for any `Atom` metadata passed in on the input system be preserved when reproducing the particles as some interatomic potential implementations may depend on this data for correctness or performance.

```julia
get_system(result::MolecularDynamicsResult, t::Integer)::AbstractSystem
```

This function should return an [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl) `AbstractSystem` which captures the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when `t` is not passed.

The default implementation creates a `FlexibleSystem` by combining the provided implementations of `get_particles`, `get_bounding_box`, and `get_boundary_conditions`.

## Simulation Analysis

```julia
temperature(result::MolecularDynamicsResult, t::Integer)::Unitful.Temperature
```

This function should return a unit-anotated temperature for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when `t` is not passed.

```julia
kinetic_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy
```

This function should return a unit-anotated kinetic energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when `t` is not passed.

```julia
potential_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy
```

This function should return a unit-anotated potential energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when `t` is not passed. If potential energy values were not cached at simulation time, the relevant function from the [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl) interface is `InteratomicPotentials.potential_energy(a::AbstractSystem, p::ArbitraryPotential)::Real`.

```julia
total_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy
```

This function should return a unit-anotated total energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when `t` is not passed.

The default implementation simply sums the kinetic and potential energy functions, but an implementor of the API might provide a custom implementation if there is a more direct means of calculation provided by the underlying simulator.

```julia
rdf(result::MolecularDynamicsResult, start::Integer, stop::Integer)::Tuple{AbstractVector{Real},AbstractVector{Real}}
```

This function should calculate the [radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function) of the system averaged across a range of timesteps. It should return a named tuple of vectors:

- r: unit-annotated interparticle radial distance bins
- g: distribution value of each bin

The default implementation provided uses the provided implementations of `get_num_bodies`, `get_bounding_box`, `get_boundary_conditions`, and `get_positions`. An implementor of the API could use a built in implementation if one that fits this spec is available.
