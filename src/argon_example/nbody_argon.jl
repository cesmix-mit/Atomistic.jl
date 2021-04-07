using LinearAlgebra
using NBodySimulator
using Plots
using StaticArrays

include("../nbody_extensions.jl")

J_per_eV = 1.602176634e-19
kg_per_amu = 1.66054e-27

kb = 1.380649e-23 # J/K

function equilibrate(N::Integer, box_size::AbstractFloat, Δt::AbstractFloat, steps::Integer, reference_temp::AbstractFloat, thermostat_prob::AbstractFloat)
    m = 6.6335209e-26 # kg
    σ = 0.34e-9 # m
    cutoff = 0.765e-9 # m
    ϵ = 1.657e-21 # J
    mean_v = √(kb * reference_temp / m) # m/s

    initial_bodies = generate_bodies_in_cell_nodes(N, m, mean_v, box_size)
	potentials = Dict(:lennard_jones => LennardJonesParameters(ϵ, σ, cutoff))
	eq_system = PotentialNBodySystem(initial_bodies, potentials)

    boundary_conditions = CubicPeriodicBoundaryConditions(box_size)
	eq_thermostat = AndersenThermostat(reference_temp, thermostat_prob / Δt)
	eq_simulation = NBodySimulation(eq_system, (0.0, steps * Δt), boundary_conditions, eq_thermostat, kb)
	
	simulator = VelocityVerlet()

    eq_result = @time run_simulation(eq_simulation, simulator, dt=Δt)
	eq_bodies = get_final_bodies(eq_result)
	return eq_result, eq_bodies
end

N = 864
# N = 100
box_size = 3.47786e-9 # m

reference_temp = 94.4 # K
thermostat_prob = 0.5

steps = 2000
Δt = 1e-14 # s

stride = 10

result, bodies = equilibrate(N, box_size, Δt, steps, reference_temp, thermostat_prob)

display(plot_temperature(result, reference_temp, stride))
display(plot_energy(result, stride))
# display(plot_rdf(result))

println(typeof(bodies))