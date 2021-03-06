# Molecular Dynamics Simulation Abstract Interface

"""
    MolecularDynamicsSimulator
Abstract type to be extended by all concrete molecular dynamics simulators.
"""
abstract type MolecularDynamicsSimulator end

"""
    simulate(system::AbstractSystem{3}, simulator::MolecularDynamicsSimulator, potential::AbstractPotential)::MolecularDynamicsResult
Run a molecular dynamics simulation configured with a particular simulator and potential with any abstract system.

An implementer of this API should implement a method of this function for their custom simulator type.
If the simulator has a fast path for some types of potential, those should be implemented with multiple dispatch.
"""
simulate(system::AbstractSystem{3}, simulator::MolecularDynamicsSimulator, potential::AbstractPotential) = throw(UnimplementedError(:simulate, simulator))
