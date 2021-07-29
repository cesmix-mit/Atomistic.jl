# Integrations with ASEPotential.jl

struct ASEPotentialParameters <: NuclearPotentialParameters
    element::ElementCoulomb
    lattice::AbstractArray{Quantity, 2}
    parameters::ASEPotential.ASECalculatorParameters
end

function forces(state::AtomicConfiguration, parameters::ASEPotentialParameters)
    atoms = ASEAtoms(state, parameters.element, parameters.lattice).atoms
    forces = get_forces(atoms, parameters.parameters)
    [@SVector [forces[i, 1], forces[i, 2], forces[i, 3]] for i ∈ 1:size(forces)[1]]
end

function potential_energy(state::AtomicConfiguration, parameters::ASEPotentialParameters)
    atoms = ASEAtoms(state, parameters.element, parameters.lattice).atoms
    get_potential_energy(atoms, parameters.parameters)
end
