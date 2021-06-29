using LinearAlgebra
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

function equilibrate(N::Integer, box_size::Quantity, Δt::Quantity, eq_steps::Integer, reference_temp::Quantity, thermostat_prob::AbstractFloat)
    m = auconvert(6.6335209e-26u"kg")
    σ = auconvert(0.34u"nm")
	cutoff = auconvert(0.765u"nm")
	ϵ = auconvert(1.657e-21u"J")
    mean_v = auconvert(√(u"k" * reference_temp / m))

    initial_bodies = generate_bodies_in_cell_nodes(N, ustrip(m), ustrip(mean_v), ustrip(box_size))
	potentials = Dict(:lennard_jones => LennardJonesParameters(ustrip(ϵ), ustrip(σ), ustrip(cutoff)))
	eq_system = PotentialNBodySystem(initial_bodies, potentials)

    boundary_conditions = CubicPeriodicBoundaryConditions(ustrip(box_size))
	eq_thermostat = AndersenThermostat(ustrip(reference_temp), thermostat_prob / ustrip(Δt))
	eq_simulation = NBodySimulation(eq_system, (0.0, eq_steps * ustrip(Δt)), boundary_conditions, eq_thermostat, 1.0)
	
	simulator = VelocityVerlet()

    return @time run_simulation(eq_simulation, simulator, dt=ustrip(Δt))
end

N = 1
box_size = auconvert(3.47786u"nm")

reference_temp = auconvert(94.4u"K")
thermostat_prob = 0.1

steps = 2000
Δt = auconvert(1e-2u"ps")

result = equilibrate(N, box_size, Δt, steps, reference_temp, thermostat_prob)
result.solution.t
