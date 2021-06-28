struct ASEPotentialParameters <: NuclearPotentialParameters
    box_size::Quantity
    element::ElementCoulomb
    lattice::AbstractArray{Quantity, 2}
    parameters::ASEPotential.ASECalculatorParameters
end

function generate_forces(bodies::AbstractVector{<:MassBody}, parameters::ASEPotentialParameters)
    atoms = dftk_atoms(parameters.element, bodies, parameters.box_size)
    atoms = ase_atoms(austrip.(parameters.lattice), atoms)
    forces = @time get_forces(atoms, parameters.parameters)
    return [@SVector [forces[i, 1], forces[i, 2], forces[i, 3]] for i ∈ 1:size(forces)[1]]
end
