using NBodySimulator
using Plots
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

function get_final_bodies(result::NBodySimulator.SimulationResult)
	N = length(result.simulation.system.bodies)
	positions = get_position(result, result.solution.t[end])
	velocities = get_velocity(result, result.solution.t[end])
	masses = get_masses(result.simulation.system)
	bodies = Array{MassBody}(undef, N)
	for i ∈ 1:N
		bodies[i] = MassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i])
	end
	return bodies
end

function plot_temperature(result::NBodySimulator.SimulationResult, stride::Integer)
    N = length(result.simulation.system.bodies)
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(result.solution.t) if (i - 1) % stride == 0]
    p = plot(
		title="Temperature during Simulation [n = $(N)]",
		xlab="Time",
		ylab="Temperature",
	)
	plot!(
		p,
		time_range,
		t -> auconvert(u"K", temperature(result, austrip(t))),
		label="Simulation Temperature",
		color=1
	)
	if (!(result.simulation.thermostat isa NBodySimulator.NullThermostat))
		plot!(
			p,
			time_range,
			t -> auconvert(u"K", result.simulation.thermostat.T),
			label="Reference Temperature",
			color=2,
			linestyle=:dash
		)
	end
	p
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
		t -> auconvert(u"hartree", kinetic_energy(result, austrip(t))),
		label="Kinetic Energy",
		color=2
	)
	plot!(
		time_range,
		t -> auconvert(u"hartree", potential_energy(result, austrip(t))),
		label="Potential Energy",
		color=1
	)
	plot!(
		time_range,
		t -> auconvert(u"hartree", total_energy(result, austrip(t))),
		label="Total Energy",
		color=3
	)
end

function plot_rdf(result::NBodySimulator.SimulationResult)
	if :lennard_jones ∉ keys(result.simulation.system.potentials)
		error("No Lennard-Jones Potential")
	end
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
		auconvert.(u"bohr", rs) / auconvert(u"bohr", σ),
		grf
	)
end

struct CustomLennardJonesParameters{pType <: Real} <: PotentialParameters
    ϵ::pType
    σ::pType
    R::pType
    σ2::pType
    R2::pType
end

function CustomLennardJonesParameters(ϵ::Real, σ::Real, R::Real)
    CustomLennardJonesParameters(ϵ, σ, R, σ^2, R^2)
end

import NBodySimulator.get_accelerating_function
function get_accelerating_function(parameters::CustomLennardJonesParameters, simulation::NBodySimulation)
    (ms, indices) = obtain_data_for_lennard_jones_interaction(simulation.system)
    (dv, u, v, t, i) -> pairwise_lennard_jones_acceleration!(dv, u, i, indices, ms, parameters, simulation.boundary_conditions)
end

function obtain_data_for_lennard_jones_interaction(system::PotentialNBodySystem)
    bodies = system.bodies
    n = length(bodies)
    ms = zeros(typeof(first(bodies).m), n)
    indices = zeros(Int, n)
    for i = 1:n
        ms[i] = bodies[i].m
        indices[i] = i
    end
    return (ms, indices)
end

function pairwise_lennard_jones_acceleration!(dv, rs, i::Integer, indices::Vector{<:Integer}, ms::Vector{<:Real}, parameters::CustomLennardJonesParameters, pbc::NBodySimulator.BoundaryConditions)
    force = @SVector [0.0, 0.0, 0.0];
    ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]

    for j ∈ indices
        if j != i
            rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            (rij, r, rij_2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

            if rij_2 < parameters.R2
                σ_rij_6 = (parameters.σ2 / rij_2)^3
                σ_rij_12 = σ_rij_6^2
                force += (2 * σ_rij_12 - σ_rij_6 ) * rij / rij_2
            end
        end
    end
    @. dv +=  24 * parameters.ϵ * force / ms[i]
end

struct DFTKForceParameters{pType <: Real} <: PotentialParameters
	forces::Array{SVector{3, pType}}
end

function get_accelerating_function(parameters::DFTKForceParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin dv .+= parameters.forces[i] / masses[i] end
end
