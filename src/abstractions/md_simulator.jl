# Molecular Dynamics Simulation Abstract Interface

abstract type MolecularDynamicsSimulator end

# TODO: state should be changed to use the common representation
function simulate(state::MassBodies, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult
    throw(UnimplementedError(:simulate, simulator))
end
