using NBodySimulator
using Plots
using StaticArrays
using UnitfulRecipes

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

function plot_temperature(result::NBodySimulator.SimulationResult, stride::Integer)
    N = length(result.simulation.system.bodies)
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(result.solution.t) if (i - 1) % stride == 0]
    plot(
		title="Temperature during Simulation [n = $(N)]",
		xlab="Time",
		ylab="Temperature",
	)
	plot!(
		time_range,
		t -> auconvert(u"K", temperature(result, ustrip(auconvert(t)))),
		label="Simulation Temperature",
		color=1
	)
	plot!(
		time_range,
		t -> auconvert(u"K", result.simulation.thermostat.T),
		label="Reference Temperature",
		color=2,
		linestyle=:dash
	)
end

function plot_energy(result::NBodySimulator.SimulationResult, stride::Integer)
    N = length(result.simulation.system.bodies)
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(result.solution.t) if (i - 1) % stride == 0]
    plot(
		title="Energy during Simulation [n = $(N)]",
		xlab="Time",
		ylab="Energy",
		legend=:right
	)
	plot!(
		time_range,
		t -> auconvert(u"hartree", kinetic_energy(result, ustrip(auconvert(t)))),
		label="Kinetic Energy",
		color=2
	)
	plot!(
		time_range,
		t -> auconvert(u"hartree", potential_energy(result, ustrip(auconvert(t)))),
		label="Potential Energy",
		color=1
	)
	plot!(
		time_range,
		t -> auconvert(u"hartree", total_energy(result, ustrip(auconvert(t)))),
		label="Total Energy",
		color=3
	)
end

function plot_rdf(result::NBodySimulator.SimulationResult)
    N = length(result.simulation.system.bodies)
    T = result.solution.destats.naccept - 1
    σ = result.simulation.system.potentials[:lennard_jones].σ
    rs, grf = @time rdf(result)
    plot(
		title="Radial Distribution Function [n = $(N)] [T = $(T)]",
		xlab="Distance r/σ",
		ylab="Radial Distribution g(r)",
        legend=false
	)
	plot!(
		[auconvert(u"bohr", r) for r ∈ rs] / auconvert(u"bohr", σ),
		grf
	)
end