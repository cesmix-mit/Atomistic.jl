module Atomistic

using DFTK
using NBodySimulator
using Plots
using PyCall
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes


include("exceptions.jl")

export MassBodies, DFTKAtoms, ASEAtoms
include("bodies.jl")

export NuclearPotentialParameters, forces, potential_energy, LJParameters
include("interatomic_potentials/abstractions.jl")
export DFTKParameters, dftk_atoms, analyze_convergence
include("interatomic_potentials/dftk_integration.jl")

export MolecularDynamicsParameters, MolecularDynamicsResult, simulate, get_bodies, get_time_range, plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf, calculate_rdf
include("simulators/abstractions.jl")
export NBSParameters, NBSResult
include("simulators/nbs_integration.jl")

export write_nbs_animation, write_ase_trajectory
include("io/visualization.jl")

end
