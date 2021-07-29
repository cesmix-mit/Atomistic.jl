# Integrations with NBodySimulator.jl

Base.@kwdef struct NBSParameters <: MolecularDynamicsParameters
    potentials::Dict{Symbol, PotentialParameters}=Dict{Symbol, PotentialParameters}()
    Δt::Quantity
    steps::Integer
    t₀::Quantity=0.0u"s"
    thermostat::NBodySimulator.Thermostat=NBodySimulator.NullThermostat()
    simulator::OrdinaryDiffEqAlgorithm=VelocityVerlet()
end

struct NBSResult <: MolecularDynamicsResult
    simulation_result::NBodySimulator.SimulationResult
end

struct CustomPotentialParameters <: PotentialParameters
    nuclear_potential_parameters::NuclearPotentialParameters
    timestep_cache::Base.RefValue{Real}
    force_cache::Base.RefValue{Vector{SVector{3, Real}}}
    CustomPotentialParameters(nuclear_potential_parameters::NuclearPotentialParameters) = new(nuclear_potential_parameters, Ref{Real}(), Ref{Vector{SVector{3, Real}}}())
end

function NBodySimulator.get_accelerating_function(parameters::CustomPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin
        if !isassigned(parameters.timestep_cache) || t != parameters.timestep_cache[]
            bodies = MassBodies(u, v, masses, simulation.boundary_conditions.L * u"bohr")
            parameters.timestep_cache[] = t
            parameters.force_cache[] = forces(bodies, parameters.nuclear_potential_parameters)
        end
        dv .+= parameters.force_cache[][i] / masses[i]
    end
end

function simulate(state::AtomicConfiguration, parameters::NBSParameters, nuclear_potential_parameters::NuclearPotentialParameters)
    parameters.potentials[:custom] = CustomPotentialParameters(nuclear_potential_parameters)
    simulate(MassBodies(state), parameters)
end

function simulate(state::AtomicConfiguration, parameters::NBSParameters, nuclear_potential_parameters::LJParameters)
    parameters.potentials[:lennard_jones] = LennardJonesParameters(austrip(nuclear_potential_parameters.ϵ), austrip(nuclear_potential_parameters.σ), austrip(nuclear_potential_parameters.R))
    simulate(MassBodies(state), parameters)
end

function simulate(state::MassBodies, parameters::NBSParameters)
    system = PotentialNBodySystem(state.bodies, parameters.potentials)
    boundary_conditions = CubicPeriodicBoundaryConditions(austrip(state.box_size))
    simulation = NBodySimulation(system, (austrip(parameters.t₀), austrip(parameters.t₀ + parameters.steps * parameters.Δt)), boundary_conditions, parameters.thermostat, 1.0)
    NBSResult(run_simulation(simulation, parameters.simulator, dt=austrip(parameters.Δt)))
end

function get_bodies(result::NBSResult, t::Integer=0)
    sr = result.simulation_result
    positions = get_position(sr, sr.solution.t[t > 0 ? t : end])
    velocities = get_velocity(sr, sr.solution.t[t > 0 ? t : end])
    masses = get_masses(sr.simulation.system)
    MassBodies(positions, velocities, masses, sr.simulation.boundary_conditions.L * u"bohr")
end

function get_time_range(result::NBSResult)
    result.simulation_result.solution.t
end

function plot_temperature!(p::Plots.Plot, result::NBSResult, stride::Integer)
    sr = result.simulation_result
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(get_time_range(result)) if (i - 1) % stride == 0]
    if (austrip(time_range[1]) != 0)
        vline!(
            p,
            [time_range[1]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        time_range,
        t -> auconvert(u"K", temperature(sr, austrip(t))),
        label=(austrip(time_range[1]) == 0 ? "Simulation Temperature" : nothing),
        color=1,
    )
    if (!(sr.simulation.thermostat isa NBodySimulator.NullThermostat))
        plot!(
            p,
            time_range,
            t -> auconvert(u"K", sr.simulation.thermostat.T),
            label=(austrip(time_range[1]) == 0 ? "Reference Temperature" : nothing),
            color=2,
            linestyle=:dash,
            lw=2
        )
    end
    p
end

function plot_energy!(p::Plots.Plot, result::NBSResult, stride::Integer)
    sr = result.simulation_result
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(get_time_range(result)) if (i - 1) % stride == 0]
    if (austrip(time_range[1]) != 0)
        vline!(
            p,
            [time_range[1]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        time_range,
        t -> kinetic_energy(sr, austrip(t))u"hartree",
        label=(austrip(time_range[1]) == 0 ? "Kinetic Energy" : nothing),
        color=2
    )
    plot!(
        p,
        time_range,
        t -> NBodySimulator.potential_energy(sr, austrip(t))u"hartree",
        label=(austrip(time_range[1]) == 0 ? "Potential Energy" : nothing),
        color=1
    )
    plot!(
        p,
        time_range,
        t -> total_energy(sr, austrip(t))u"hartree",
        label=(austrip(time_range[1]) == 0 ? "Total Energy" : nothing),
        color=3
    )
end

function calculate_rdf(result::NBSResult, sample_fraction::Real)
    sr = result.simulation_result
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    trange = get_time_range(result)[end - floor(Int, length(sr.solution.t) * sample_fraction) + 1:end]

    maxbin = 1000
    dr = pbc.L / 2 / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for i ∈ 1:n
            ri = @SVector [cc[1, i], cc[2, i], cc[3, i]]
            for j ∈ i + 1:n
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

    c = 4 / 3 * π * n / pbc.L^3
    gr = zeros(maxbin)
    rs = zeros(maxbin)
    tlen = length(trange)
    for bin ∈ 1:maxbin
        rlower = (bin - 1) * dr
        rupper = rlower + dr
        nideal = c * (rupper^3 - rlower^3)
        gr[bin] = (hist[bin] / (tlen * n)) / nideal
        rs[bin] = rlower + dr / 2
    end

    rs, gr
end
