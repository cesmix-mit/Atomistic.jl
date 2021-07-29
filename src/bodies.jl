# Body Collection Abstractions

abstract type AtomicConfiguration end

struct MassBodies <: AtomicConfiguration
    bodies::Vector{MassBody}
    box_size::Quantity
end

struct DFTKAtoms <: AtomicConfiguration
    atoms
    box_size::Quantity
    lattice::AbstractArray{Quantity, 2}
end

struct ASEAtoms <: AtomicConfiguration
    atoms
    box_size::Quantity
end

MassBodies(state::MassBodies) = state
function MassBodies(positions::AbstractMatrix{<:Real}, velocities::AbstractMatrix{<:Real}, masses::AbstractVector{<:Real}, box_size::Quantity)
    positions = mod.(positions, austrip(box_size)) # https://github.com/SciML/NBodySimulator.jl/issues/46
    MassBodies([MassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i]) for i ∈ 1:length(masses)], box_size)
end
function MassBodies(N::Integer, mass::Quantity, average_velocity::Quantity, box_size::Quantity)
    MassBodies(generate_bodies_in_cell_nodes(N, austrip(mass), austrip(average_velocity), austrip(box_size)), box_size)
end

DFTKAtoms(state::DFTKAtoms) = state
DFTKAtoms(state::ASEAtoms) = DFTKAtoms(load_atoms_ase(state.atoms), state.box_size, load_lattice(state.atoms)u"bohr")
function DFTKAtoms(state::MassBodies, element::DFTK.Element, lattice::AbstractArray{<:Quantity, 2})
    DFTKAtoms([element => [austrip.(b.r * u"bohr" / state.box_size) for b ∈ state.bodies]], state.box_size, lattice)
end

ASEAtoms(state::ASEAtoms) = state
ASEAtoms(state::DFTKAtoms) = ASEAtoms(ase_atoms(austrip.(state.lattice), state.atoms), state.box_size)
function ASEAtoms(state::MassBodies, element::DFTK.Element, lattice::AbstractArray{<:Quantity, 2})
    ASEAtoms(DFTKAtoms(state, element, lattice))
end

Base.length(state::MassBodies) = length(state.bodies)
Base.length(state::DFTKAtoms) = length(state.atoms)
