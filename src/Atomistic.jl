module Atomistic

using AtomsBase
using InteratomicPotentials
using DFTK
using LinearAlgebra
using NBodySimulator
using PeriodicTable
using Plots
using PyCall
using Random
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

import Base: @kwdef, Fix2
import InteratomicPotentials: ArbitraryPotential
import DFTK: Mixing
import NBodySimulator: Body, BoundaryConditions, NullThermostat, SimulationResult, Thermostat
import Plots: Plot

include("unit_convention.jl")

# AtomsBase Integrations
export DynamicSystem
include("integrations/atomsbase_integration.jl")

# DFTK Integrations
export DFTKPotential
export analyze_convergence
include("integrations/dftk_integration.jl")

# API
include("api/exceptions.jl")
export MolecularDynamicsResult
export get_time_range, get_bounding_box, get_boundary_conditions, reference_temperature
export get_time, get_positions, get_velocities, get_particles, get_system
export temperature, kinetic_energy, potential_energy, total_energy, rdf
include("api/md_result.jl")
export MolecularDynamicsSimulator
export simulate
include("api/md_simulator.jl")

# Implementations
export NBSimulator, NBSResult
include("implementations/nbodysimulator.jl")

# Analysis
export plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf
include("analysis/plotting.jl")
export write_nbs_animation, write_ase_trajectory
include("analysis/visualization.jl")

end
