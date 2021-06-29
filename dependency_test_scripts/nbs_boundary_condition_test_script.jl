using NBodySimulator
using Unitful
using UnitfulAtomic

m = 6.6335209e-26u"kg"
σ = 0.34u"nm"
cutoff = 0.765u"nm"
ϵ = 1.657e-21u"J"

argon_initial_bodies(N::Integer, box_size::Quantity, reference_temp::Quantity) = generate_bodies_in_cell_nodes(N, austrip(m), austrip(√(u"k" * reference_temp / m)), austrip(box_size))
argon_lennard_jones() = Dict(:lennard_jones => LennardJonesParameters(austrip(ϵ), austrip(σ), austrip(cutoff)))
argon_equilibration_thermostat(reference_temp::Quantity, thermostat_prob::AbstractFloat, Δt::Quantity) = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))

N = 8
box_size = 20u"bohr"
reference_temp = 94.4u"K"
thermostat_prob = 0.1
t₀ = 0u"ps"
Δt = 1e-2u"ps"
steps = 20000

initial_bodies = argon_initial_bodies(N, box_size, reference_temp)
potentials = argon_lennard_jones()
system = PotentialNBodySystem(initial_bodies, potentials)
boundary_conditions = CubicPeriodicBoundaryConditions(austrip(box_size))
thermostat = argon_equilibration_thermostat(reference_temp, thermostat_prob, Δt)
simulation = NBodySimulation(system, (austrip(t₀), austrip(t₀ + steps * Δt)), boundary_conditions, thermostat, 1.0)
simulator = VelocityVerlet()

result = run_simulation(simulation, simulator, dt=austrip(Δt))

positions = get_position(result, result.solution.t[end])

if isa(boundary_conditions, CubicPeriodicBoundaryConditions)
    bounded_positions = mod.(positions, boundary_conditions.L)
end

display(positions)
display(bounded_positions)
