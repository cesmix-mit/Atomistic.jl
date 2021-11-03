# Molecular Dynamics Simulation Abstract Interface

abstract type MolecularDynamicsSimulator end

function simulate(system::AbstractSystem, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult
    throw(UnimplementedError(:simulate, simulator))
end
