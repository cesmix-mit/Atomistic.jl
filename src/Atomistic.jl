module Atomistic

using AtomsBase
using InteratomicPotentials
using DFTK
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

# API
include("api/exceptions.jl")
export MolecularDynamicsResult, get_system, get_time_range, temperature, reference_temperature, kinetic_energy, potential_energy, total_energy, rdf
include("api/md_result.jl")
export MolecularDynamicsSimulator, simulate
include("api/md_simulator.jl")

# Integrations
export DynamicSystem
include("integrations/atomsbase_integration.jl")
export DFTKPotential, analyze_convergence
include("integrations/dftk_integration.jl")

# Implementations
export NBSimulator, NBSResult
include("implementations/nbodysimulator.jl")

# Analysis
export plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf
include("analysis/plotting.jl")
export write_nbs_animation, write_ase_trajectory
include("analysis/visualization.jl")

end
