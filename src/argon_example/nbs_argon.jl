using LinearAlgebra
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

include("../nbs_extensions.jl")

function equilibrate(N::Integer, box_size::Quantity, Δt::Quantity, eq_steps::Integer, reference_temp::Quantity, thermostat_prob::AbstractFloat)
    m = auconvert(6.6335209e-26u"kg")
    σ = auconvert(0.34u"nm")
	cutoff = auconvert(0.765u"nm")
	ϵ = auconvert(1.657e-21u"J")
    mean_v = auconvert(√(u"k" * reference_temp / m))

    initial_bodies = generate_bodies_in_cell_nodes(N, ustrip(m), ustrip(mean_v), ustrip(box_size))
	potentials = Dict(:lennard_jones => LennardJonesParameters(ustrip(ϵ), ustrip(σ), ustrip(cutoff)))
	# potentials = Dict(:custom_lennard_jones => CustomLennardJonesParameters(ustrip(ϵ), ustrip(σ), ustrip(cutoff)))
	eq_system = PotentialNBodySystem(initial_bodies, potentials)

    boundary_conditions = CubicPeriodicBoundaryConditions(ustrip(box_size))
	eq_thermostat = AndersenThermostat(ustrip(reference_temp), thermostat_prob / ustrip(Δt))
	eq_simulation = NBodySimulation(eq_system, (0.0, eq_steps * ustrip(Δt)), boundary_conditions, eq_thermostat, 1.0)
	
	simulator = VelocityVerlet()

    eq_result = @time run_simulation(eq_simulation, simulator, dt=ustrip(Δt))
	eq_bodies = get_final_bodies(eq_result)
	return eq_result, eq_bodies
end

function dftk_force_simulate(bodies::Vector{MassBody}, forces::Vector, box_size::Quantity, Δt::Quantity, steps::Integer)
	potentials = Dict(:DFTK => DFTKForceParameters([ustrip.(f) for f in forces]))
	system = PotentialNBodySystem(bodies, potentials)

	boundary_conditions = CubicPeriodicBoundaryConditions(ustrip(box_size))
	simulation = NBodySimulation(system, (0.0, steps * ustrip(Δt)), boundary_conditions, 1.0)

	simulator = VelocityVerlet()

	result = @time run_simulation(simulation, simulator, dt=ustrip(Δt))
	bodies = get_final_bodies(result)
	return result, bodies
end
