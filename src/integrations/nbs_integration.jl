# Integrations with NBodySimulator.jl

@kwdef struct NBSimulator <: MolecularDynamicsSimulator
    Δt::Real
    steps::Integer
    t₀::Real = 0.0
    thermostat::Thermostat = NullThermostat()
    simulator::OrdinaryDiffEqAlgorithm = VelocityVerlet()
    potentials::Dict{Symbol,PotentialParameters} = Dict{Symbol,PotentialParameters}()
end
function NBSimulator(;
                     Δt::Quantity,
                     steps::Integer,
                     t₀::Quantity=0.0u"s",
                     thermostat::Thermostat=NullThermostat(),
                     simulator::OrdinaryDiffEqAlgorithm=VelocityVerlet(),
                     potentials::Dict{Symbol,PotentialParameters}=Dict{Symbol,PotentialParameters}())
    NBSimulator(austrip(Δt), steps, austrip(t₀), thermostat, simulator, potentials)
end

@kwdef struct CustomPotentialParameters <: PotentialParameters
    potential::ArbitraryPotential
    timestep_cache::RefValue{Real} = Ref{Real}()
    force_cache::RefValue{Vector{SVector{3,Real}}} = Ref{Vector{SVector{3,Real}}}()
end

@kwdef struct LJPotential
    ϵ::Real
    σ::Real
    R::Real
end
function LJPotential(; ϵ::Quantity, σ::Quantity, R::Quantity)
    LJPotential(austrip(ϵ), austrip(σ), austrip(R))
end

function NBodySimulator.get_accelerating_function(parameters::CustomPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin
        if !isassigned(parameters.timestep_cache) || t != parameters.timestep_cache[]
            bodies = MassBodies(u, v, masses, simulation.boundary_conditions.L * u"bohr")
            parameters.timestep_cache[] = t
            parameters.force_cache[] = force(bodies, parameters.potential)
        end
        dv .+= parameters.force_cache[][i] / masses[i]
    end
end

function simulate(state::MassBodies, simulator::NBSimulator, potential::ArbitraryPotential)
    simulator.potentials[:custom] = CustomPotentialParameters(potential=potential)
    simulate(state, simulator)
end

function simulate(state::MassBodies, simulator::NBSimulator, potential::LJPotential)
    simulator.potentials[:lennard_jones] = LennardJonesParameters(potential.ϵ, potential.σ, potential.R)
    simulate(state, simulator)
end

function simulate(state::MassBodies, simulator::NBSimulator)
    system = PotentialNBodySystem(state.bodies, simulator.potentials)
    boundary_conditions = CubicPeriodicBoundaryConditions(austrip(state.box_size))
    simulation = NBodySimulation(system, (simulator.t₀, simulator.t₀ + simulator.steps * simulator.Δt), boundary_conditions, simulator.thermostat, 1.0)
    NBSResult(run_simulation(simulation, simulator.simulator, dt=simulator.Δt))
end

struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
end

function get_bodies(result::NBSResult, t::Integer=0)
    sr = result.result
    positions = get_position(sr, sr.solution.t[t > 0 ? t : end])
    velocities = get_velocity(sr, sr.solution.t[t > 0 ? t : end])
    masses = get_masses(sr.simulation.system)
    MassBodies(positions, velocities, masses, sr.simulation.boundary_conditions.L * u"bohr")
end

function get_time_range(result::NBSResult)
    result.result.solution.t
end

function temperature(result::NBSResult, time::Real)
    NBodySimulator.temperature(result.result, time)
end

function reference_temperature(result::NBSResult)
    result.result.simulation.thermostat isa NullThermostat ? nothing : result.result.simulation.thermostat.T
end

function kinetic_energy(result::NBSResult, time::Real)
    NBodySimulator.kinetic_energy(result.result, time)
end

function potential_energy(result::NBSResult, time::Real)
    NBodySimulator.potential_energy(result.result, time)
end

function total_energy(result::NBSResult, time::Real)
    NBodySimulator.total_energy(result.result, time)
end

function rdf(result::NBSResult, sample_fraction::Float64=1.0)
    @assert 0 < sample_fraction ≤ 1
    sr = result.result
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
