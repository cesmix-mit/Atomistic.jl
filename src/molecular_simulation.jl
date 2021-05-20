using DFTK
using NBodySimulator
using Plots
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

abstract type ForceGenerationParameters end

struct ParticleForcePotentialParameters <: PotentialParameters
	forces::Vector{SVector{3, Real}}
end
ParticleForcePotentialParameters(forces::Matrix{<:Real}) = ParticleForcePotentialParameters([@SVector [forces[i, 1], forces[i, 2], forces[i, 3]] for i ∈ 1:size(forces)[1]])

function NBodySimulator.get_accelerating_function(parameters::ParticleForcePotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin dv .+= parameters.forces[i] / masses[i] end
end
