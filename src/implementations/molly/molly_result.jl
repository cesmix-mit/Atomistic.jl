# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsResult API
# -----------------------------------------------------------------------------

"""
    MollyResult <: MolecularDynamicsResult

The result generated from running an `MollySimulator`.

**Field descriptions**
- `system::System` the `System` resulting from the simulation
- `simulator::MollySimulator` the `MollySimulator` that produced the result
"""
struct MollyResult <: MolecularDynamicsResult
    system::System
    simulator::MollySimulator
end

function get_time_range(result::MollyResult)
    s = result.simulator
    Δt = s.Δt * s.stride
    start = s.t₀ + Δt - s.Δt
    steps = s.steps ÷ s.stride - 1
    start:Δt:start+Δt*steps
end
get_num_bodies(result::MollyResult) = length(result.system)
get_bounding_box(result::MollyResult) = bounding_box(result.system)
get_boundary_conditions(result::MollyResult) = boundary_conditions(result.system)

reference_temperature(result::MollyResult) = coupling_reference_temperature(result.simulator.coupling)
coupling_reference_temperature(coupling) = coupling.temperature
coupling_reference_temperature(::NoCoupling) = missing

get_positions(result::MollyResult) = position(result.system)
get_velocities(result::MollyResult) = velocity(result.system)
get_particles(result::MollyResult) = AtomsBase.Atom.(result.system.atoms_data, get_positions(result), get_velocities(result))
get_system(result::MollyResult) = result.system

get_positions(result::MollyResult, t::Integer) = result.system.loggers["c"].coords[t]
get_velocities(result::MollyResult, t::Integer) = result.system.loggers["v"].velocities[t]
get_particles(result::MollyResult, t::Integer) = AtomsBase.Atom.(result.system.atoms_data, get_positions(result, t), get_velocities(result, t))

temperature(result::MollyResult, t::Integer) = result.system.loggers["t"].temperatures[t]
kinetic_energy(result::MollyResult, t::Integer) = result.system.loggers["k"].energies[t]
potential_energy(result::MollyResult, t::Integer) = result.system.loggers["p"].energies[t]