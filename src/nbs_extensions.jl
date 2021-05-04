using NBodySimulator
using Plots
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes

function simulate(bodies::Vector{MassBody}, potential::PotentialParameters, box_size::Quantity, Δt::Quantity, steps::Integer)
	potentials = Dict(:custom => potential)
	system = PotentialNBodySystem(bodies, potentials)

	boundary_conditions = CubicPeriodicBoundaryConditions(ustrip(box_size))
	simulation = NBodySimulation(system, (0.0, steps * ustrip(Δt)), boundary_conditions, 1.0)

	simulator = VelocityVerlet()

	result = @time run_simulation(simulation, simulator, dt=ustrip(Δt))
	bodies = get_final_bodies(result)
	return result, bodies
end

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
    T = result.solution.destats.naccept
	σ = result.simulation.system.potentials[:lennard_jones].σ
    rs, grf = @time custom_rdf(result, 20)
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

import NBodySimulator.obtain_data_for_lennard_jones_interaction
function custom_rdf(sr::NBodySimulator.SimulationResult, sample_fraction::Integer=10)
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions

    (ms, indxs) = obtain_data_for_lennard_jones_interaction(sr.simulation.system)
    indlen = length(indxs)
	trange = sr.solution.t[end - length(sr.solution.t) ÷ sample_fraction + 1:end]

    maxbin = 1000
    dr = pbc.L / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for ind_i = 1:indlen
            i = indxs[ind_i]
            ri = @SVector [cc[1, i], cc[2, i], cc[3, i]]
            for ind_j = ind_i + 1:indlen
                j = indxs[ind_j]
                rj = @SVector [cc[1, j], cc[2, j], cc[3, j]]

                (rij, r, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)

                if r2 < (0.5 * pbc.L)^2
                    bin = ceil(Int, r / dr)
                    if bin > 1 && bin <= maxbin
                        hist[bin] += 2
                    end
                end
            end
        end
    end

    c = 4 / 3 * π * indlen / pbc.L^3

    gr = zeros(maxbin)
    rs = zeros(maxbin)
    tlen = length(trange)
    for bin = 1:maxbin
        rlower = (bin - 1) * dr
        rupper = rlower + dr
        nideal = c * (rupper^3 - rlower^3)
        gr[bin] = (hist[bin] / (tlen * indlen)) / nideal
        rs[bin] = rlower + dr / 2
    end

    return (rs, gr)
end
