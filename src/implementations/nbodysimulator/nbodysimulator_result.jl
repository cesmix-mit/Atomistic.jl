# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsResult API
# -----------------------------------------------------------------------------

"""
    NBSResult <: MolecularDynamicsResult

The result generated from running an `NBSimulator`.

**Field descriptions**
- `result::SimulationResult` the standard simulation result from `NBodySimulator`
- `energy_cache::Vector{Float64}` cached potential energy values from the simulation
"""
struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
    energy_cache::Vector{Float64}
end

get_time_range(result::NBSResult) = result.result.solution.t * TIME_UNIT
get_num_bodies(result::NBSResult) = length(result.result.simulation.system.bodies)
get_bounding_box(result::NBSResult) = get_bounding_box(result.result.simulation.boundary_conditions)
get_boundary_conditions(result::NBSResult) = get_boundary_conditions(result.result.simulation.boundary_conditions)

reference_temperature(result::NBSResult) = reference_temperature(result.result.simulation.thermostat)
reference_temperature(thermostat::Thermostat) = thermostat.T * TEMPERATURE_UNIT
reference_temperature(::NullThermostat) = missing

function get_positions(result::NBSResult, t::Integer)
    positions = get_position(result.result, result.result.solution.t[t])
    [bound_position(SVector{3}(p), result.result.simulation.boundary_conditions) * LENGTH_UNIT for p ∈ eachcol(positions)]
end
function get_velocities(result::NBSResult, t::Integer)
    velocities = get_velocity(result.result, result.result.solution.t[t])
    [SVector{3}(v) * VELOCITY_UNIT for v ∈ eachcol(velocities)]
end
get_particles(result::NBSResult, t::Integer) = AtomsBase.Atom.(result.result.simulation.system.bodies, get_positions(result, t), get_velocities(result, t))

temperature(result::NBSResult, t::Integer) = NBodySimulator.temperature(result.result, result.result.solution.t[t]) * TEMPERATURE_UNIT
kinetic_energy(result::NBSResult, t::Integer) = NBodySimulator.kinetic_energy(result.result, result.result.solution.t[t]) * ENERGY_UNIT
potential_energy(result::NBSResult, t::Integer) = result.energy_cache[t] * ENERGY_UNIT

"""
    animate(result::NBSResult, filename::String; kwargs...)

Animate an `NBSResult` using the interface provided by NBodySimulator.jl and store the result in a file.
"""
Plots.animate(result::NBSResult, filename::String; kwargs...) = animate(result.result, filename; kwargs...)
