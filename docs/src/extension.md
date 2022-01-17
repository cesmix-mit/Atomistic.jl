# Implementing the Atomistic API

There are two main parts to the Atomistic API. An implementer must provide a concrete simulator type subtying the abstract type `MolecularDynamicsSimulator` and a concrete result type subtyping the abstract type `MolecularDynamicsResult`. A summary of the required functionality for each component follows.

# Molecular Dynamics Simulators

The only function that must be implemented for subtypes of `MolecularDynamicsSimulator` is a method of `simulate`. The signature for this function is the following

```julia
simulate(system::AbstractSystem, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult
```

The first parameter is an instance of an implementation of the `AbstractSystem` interface defined by [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl). It is important that the implementer only rely on this API rather than any particular implementation details because the user could provide any system implementation. The species, position, and velocity of each body as well as the boundary shape and conditions are stored in this data structure. An implementer should take special care to handle units appropriately when converting from this representation to the internal representation of the underlying simulator package.

The second parameter is the user-defined simulator. This could specify information such as the duration of the simulation, numerical details such as discretization parameters and integration methods, or domain specific options such as thermostats and barostats.

The third parameter is an instance of an implementation of the `ArbitraryPotential` interface defined by [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl). The main function from this interface that is useful in simulation context is `force(s::AbstractSystem, p::ArbitraryPotential)::AbstractVector{StaticVector{3, <:Real}}`.

The `simulate` method should return the corresponding implementation of the `MolecularDynamicsResult` interface as described below.

# [Molecular Dynamics Results](@id MolecularDynamicsResult_Specification)

Eight total functions must be implemented for subtypes of `MolecularDynamicsResult`, two of which have default implementations. The functions allow users to access simulation data from three categories: simulation configuration, time-series measurable quantities, and simulation analysis.

## Simulation Configuration Functions

There are two simulation configuration function for `MolecularDynamicsResult`s:

```julia
get_time_range(result::MolecularDynamicsResult)::AbstractVector{<:Unitful.Time}
```

This function should return a unit-anotated vector-like object containing the time value for each step of the simulation which can be iterated for plotting, animation, or other analysis.

```julia
reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}
```

This function should return a unit-anotated temperature that describes the thermostat used in the simulation. If there is no reference temperature for the particular simulation, the function should return missing, which is the default implementation.

## Time-Series Measureable Quantities

There are 5 functions in the `MolecularDynamicsResult` API which return measured quantities from a particular timestep.

```julia
get_system(result::MolecularDynamicsResult, t::Integer = 0)::AbstractSystem
```

This function should return an [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl) `AbstractSystem` which describes the system at a particular timestep in the simulation. An implementer should take special care that any unit transformations in this stage are done appropriately.

```julia
temperature(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Temperature
```

This function should return a unit-anotated temperature for the system at a particular timestep in the simulation.

```julia
kinetic_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy
```

This function should return a unit-anotated kinetic energy for the system at a particular timestep in the simulation.

```julia
potential_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy
```

This function should return a unit-anotated potential energy for the system at a particular timestep in the simulation. The relevant function from the [InteratomicPotentials.jl](<(https://github.com/cesmix-mit/InteratomicPotentials.jl)>) interface is `potential_energy(a::AbstractSystem, p::ArbitraryPotential)::Real`.

```julia
total_energy(result::MolecularDynamicsResult, t::Integer = 0)::Unitful.Energy
```

This function should return a unit-anotated total energy for the system at a particular timestep in the simulation. The default implementation simply sums the kinetic and potential energy functions, but an implemention might provide a custom implementation if there is a more direct means of calculation provided by the underlying simulator.

## Simulation Analysis

The only analysis function required by the Atmostic interface is the [Radial Distribution Function](https://en.wikipedia.org/wiki/Radial_distribution_function).

```julia
rdf(result::MolecularDynamicsResult, sample_fraction::Float64 = 1.0)::Tuple{AbstractVector{<:Real},AbstractVector{<:Real}}
```

This function should calculate a tuple of vectors which represent the interparticle radial distances (in bohr) and the density of each distance respectively. The densities should be averaged across a trailing fraction of timesteps. Note that this detail is still the subject of further scrutiny and might be modified in a future release.
