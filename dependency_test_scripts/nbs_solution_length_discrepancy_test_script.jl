using NBodySimulator
using Unitful
using UnitfulAtomic

m = 6.6335209e-26u"kg"
σ = 0.34u"nm"
cutoff = 0.765u"nm"
ϵ = 1.657e-21u"J"

argon_initial_bodies(N::Integer, box_size::Quantity, reference_temp::Quantity) = generate_bodies_in_cell_nodes(N, austrip(m), austrip(√(u"k" * reference_temp / m)), austrip(box_size))
argon_lennard_jones() = Dict(:lennard_jones => LennardJonesParameters(austrip(ϵ), austrip(σ), austrip(cutoff)))
argon_equilibration_thermostat(reference_temp::Quantity, thermostat_prob::AbstractFloat, Δt::AbstractFloat) = AndersenThermostat(austrip(reference_temp), thermostat_prob / Δt)

N = 8
box_size = 20u"bohr"
reference_temp = 94.4u"K"
thermostat_prob = 0.1
t₀ = 0.0
Δt = 413.41373335335163

initial_bodies = argon_initial_bodies(N, box_size, reference_temp)
potentials = argon_lennard_jones()
system = PotentialNBodySystem(initial_bodies, potentials)
boundary_conditions = CubicPeriodicBoundaryConditions(austrip(box_size))
thermostat = argon_equilibration_thermostat(reference_temp, thermostat_prob, Δt)
simulator = VelocityVerlet()

for steps ∈ (10000, 20000, 30000)
    simulation = NBodySimulation(system, (t₀, t₀ + steps * Δt), boundary_conditions, thermostat, 1.0)

    result = run_simulation(simulation, simulator, dt=Δt)

    display(length(result.solution.t))
    # expected output: 10001, 20001, 30001
    # actual output: 10001, 20002, 30001
end
