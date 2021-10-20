module Atomistic

using InteratomicPotentials
using DFTK
using NBodySimulator
using Plots
using PyCall
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

import Base:@kwdef, RefValue
import InteratomicPotentials:ArbitraryPotential
import DFTK:Mixing, Element
import NBodySimulator:Thermostat, NullThermostat, SimulationResult
import Plots:Plot

include("exceptions.jl")

export MassBodies, DFTKAtoms, ASEAtoms
include("bodies.jl")

export MolecularDynamicsResult, get_bodies, get_time_range, rdf
include("abstractions/md_result.jl")
export MolecularDynamicsSimulator, simulate
include("abstractions/md_simulator.jl")

include("integrations/ip_integration.jl")
export NBSimulator, LJParameters, NBSResult
include("integrations/nbs_integration.jl")
export DFTKPotential, analyze_convergence
include("integrations/dftk_integration.jl")

export plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf
include("analysis/ploting.jl")
export write_nbs_animation, write_ase_trajectory
include("analysis/visualization.jl")

end
