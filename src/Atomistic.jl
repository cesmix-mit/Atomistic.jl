module Atomistic

using Reexport

using AtomsBase
using InteratomicPotentials
using Distances
using LinearAlgebra
using StaticArrays
@reexport using Unitful
@reexport using UnitfulAtomic

using DFTK
using Molly
using NBodySimulator

using Random
using PeriodicTable
using Plots
using PyCall
using UnitfulRecipes

import Base: @kwdef
import InteratomicPotentials: ArbitraryPotential
import DFTK: Mixing
import Molly: AbstractNeighborFinder
import NBodySimulator: Body, BoundaryConditions, NullThermostat, SimulationResult, Thermostat
import Plots: Plot

include("unit_convention.jl")

# API
include("api/exceptions.jl")
export MolecularDynamicsResult
export get_time_range, get_num_bodies, get_bounding_box, get_boundary_conditions, reference_temperature
export get_time, get_positions, get_velocities, get_particles, get_system
export temperature, kinetic_energy, potential_energy, total_energy, rdf
include("api/md_result.jl")
export MolecularDynamicsSimulator
export simulate
include("api/md_simulator.jl")

# Implementations
export NBSimulator, NBSResult
include("implementations/nbodysimulator.jl")
export MollySimulator, MollyResult
include("implementations/molly.jl")

# Initialization Convenience functions
export generate_atoms_in_cubic_cell
include("integrations/initialization.jl")

# Analysis
export plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf
include("analysis/plotting.jl")
export write_nbs_animation, write_ase_trajectory
include("analysis/visualization.jl")

# DFTK Integrations
export DFTKPotential
include("integrations/dftk_integration.jl")

end
