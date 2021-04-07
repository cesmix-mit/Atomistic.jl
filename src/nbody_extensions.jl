using NBodySimulator
using Plots
using StaticArrays

function get_final_bodies(result::NBodySimulator.SimulationResult)
	N = length(result.simulation.system.bodies)
	positions = get_position(result, result.solution.t[end])
	velocities = get_velocity(result, result.solution.t[end])
	masses = get_masses(result.simulation.system)
	bodies = MassBody[]
	for i ∈ 1:N
		push!(bodies, MassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i]))
	end
	return bodies
end

function plot_temperature(result::NBodySimulator.SimulationResult, reference_temp::AbstractFloat, stride::Integer, time_scale::AbstractFloat=1e12)
    N = length(result.simulation.system.bodies)
    reference_temp = result.simulation.thermostat.T
    time_range = [t * time_scale for (i, t) ∈ enumerate(result.solution.t) if (i - 1) % stride == 0]
    plot(
		title="Temperature during Simulation [n = $(N)]",
		xlab="Time [ps]",
		ylab="Temperature [K]",
	)
	plot!(
		time_range,
		t -> temperature(result, t / time_scale),
		label="Simulation Temperature",
		color=1
	)
	plot!(
		time_range,
		t -> reference_temp,
		label="Reference Temperature",
		color=2,
		linestyle=:dash
	)
end

function plot_energy(result::NBodySimulator.SimulationResult, stride::Integer, time_scale::AbstractFloat=1e12)
    N = length(result.simulation.system.bodies)
    time_range = [t * time_scale for (i, t) ∈ enumerate(result.solution.t) if (i - 1) % stride == 0]
    plot(
		title="Energy during Simulation [n = $(N)]",
		xlab="Time [ps]",
		ylab="Energy [eV]",
		legend=:right
	)
	plot!(
		time_range,
		t -> kinetic_energy(result, t * 1e-12) / J_per_eV,
		label="Kinetic Energy",
		color=2
	)
	plot!(
		time_range,
		t -> potential_energy(result, t * 1e-12) / J_per_eV,
		label="Potential Energy",
		color=1
	)
	plot!(
		time_range,
		t -> total_energy(result, t * 1e-12) / J_per_eV,
		label="Total Energy",
		color=3
	)
end

function plot_rdf(result::NBodySimulator.SimulationResult)
    N = length(result.simulation.system.bodies)
    σ = result.simulation.system.potentials[:lennard_jones].σ
    rs, grf = @time rdf(result)
    plot(
		title="Radial Distribution Function [n = $(N)]",
		xlab="Distance r/σ",
		ylab="Radial Distribution g(r)",
        legend=false
	)
	plot!(
		rs / σ,
		grf
	)
end