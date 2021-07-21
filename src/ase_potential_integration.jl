# Integrations with ASEPotential.jl

struct ASEPotentialParameters <: NuclearPotentialParameters
    box_size::Quantity
    element::ElementCoulomb
    lattice::AbstractArray{Quantity, 2}
    parameters::ASEPotential.ASECalculatorParameters
end

function generate_forces(bodies::AbstractVector{<:MassBody}, parameters::ASEPotentialParameters)
    atoms = ase_atoms(parameters.element, bodies, parameters.box_size, parameters.lattice)
    forces = @time get_forces(atoms, parameters.parameters)
    return [@SVector [forces[i, 1], forces[i, 2], forces[i, 3]] for i ∈ 1:size(forces)[1]]
end

function ase_atoms(element::DFTK.Element, bodies::AbstractVector{<:MassBody}, box_size::Quantity, lattice::AbstractArray{<:Quantity, 2})
    DFTK.ase_atoms(austrip.(lattice), dftk_atoms(element, bodies, box_size))
end
