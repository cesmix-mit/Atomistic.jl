# Integrations with ASEPotential.jl

struct ASEPotentialParameters <: NuclearPotentialParameters
    element::ElementCoulomb
    lattice::AbstractArray{Quantity, 2}
    parameters::ASEPotential.ASECalculatorParameters
end

function generate_forces(state::AtomicConfiguration, parameters::ASEPotentialParameters)
    atoms = ASEAtoms(state, parameters.element, parameters.lattice).atoms
    forces = get_forces(atoms, parameters.parameters)
    return [@SVector [forces[i, 1], forces[i, 2], forces[i, 3]] for i ∈ 1:size(forces)[1]]
end
