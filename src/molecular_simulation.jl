using ASEPotential
using DFTK
using NBodySimulator
using Plots
using PyCall
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

include("nbs_extensions.jl")
include("dftk_integration.jl")
include("ase_potential_integration.jl")
include("io/ase_trajectory.jl")

abstract type ForceGenerationParameters <: PotentialParameters end

function NBodySimulator.get_accelerating_function(parameters::ForceGenerationParameters, simulation::NBodySimulation)
    forces = generate_forces(simulation.system.bodies, parameters)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin dv .+= forces[i] / masses[i] end
end

Base.@kwdef struct AbInitioPotentialParameters <: PotentialParameters
    forceGenerationParameters::ForceGenerationParameters
    forceCache::Dict{Real, Vector{SVector{3, Real}}} = Dict{Real, Vector{SVector{3, Real}}}()
end

function NBodySimulator.get_accelerating_function(parameters::AbInitioPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin
        if t ∉ keys(parameters.forceCache)
            bodies = construct_bodies(u, v, masses, simulation.boundary_conditions)
            parameters.forceCache[t] = generate_forces(bodies, parameters.forceGenerationParameters)
        end
        dv .+= parameters.forceCache[t][i] / masses[i]
    end
end
