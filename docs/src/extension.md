# Implementing the Atomistic API

TODO

## [API Reference](@id API_Reference)

### Molecular Dynamics Simulator

```@docs
MolecularDynamicsSimulator
simulate(::AbstractSystem, ::MolecularDynamicsSimulator, ::ArbitraryPotential)
```

### Molecular Dynamics Result

```@docs
MolecularDynamicsResult
get_system(::MolecularDynamicsResult, ::Integer)
get_time_range(::MolecularDynamicsResult)
temperature(::MolecularDynamicsResult, ::Integer)
reference_temperature(::MolecularDynamicsResult)
kinetic_energy(::MolecularDynamicsResult, ::Integer)
potential_energy(::MolecularDynamicsResult, ::Integer)
total_energy(::MolecularDynamicsResult, ::Integer)
rdf(::MolecularDynamicsResult, ::Float64)
```
