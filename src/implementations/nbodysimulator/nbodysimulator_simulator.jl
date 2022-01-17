# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsSimulator API
# -----------------------------------------------------------------------------

"""
    NBSimulator <: MolecularDynamicsSimulator

A wrapper around NBodySimulator to implement the Atomistic API.

**Field descriptions**
- `Δt::Real` the time between timesteps, assumed to be in atomic units
- `steps::Integer` the number of timesteps for the simulation
- `t₀::Real` the starting time of the simulation, assumed to be in atomic units;
    defaults to 0
- `thermostat::Thermostat` the thermostat for the simulation;
    many options are defined by `NBodySimulator`, but a user could also define a custom thermostat;
    defaults to the `NullThermostat`
- `simulator::OrdinaryDiffEqAlgorithm` the algorithm to be used for the ODE;
    defaults to VelocityVerlet
- `potentials::Dict{Symbol,PotentialParameters}` dictionary of potentials;
    shouldn't be manipulated directly by the user
"""
@kwdef struct NBSimulator <: MolecularDynamicsSimulator
    Δt::Real       # in TIME_UNIT
    steps::Integer
    t₀::Real = 0.0 # in TIME_UNIT
    thermostat::Thermostat = NullThermostat()
    simulator::OrdinaryDiffEqAlgorithm = VelocityVerlet()
    potentials::Dict{Symbol,PotentialParameters} = Dict{Symbol,PotentialParameters}()
end
function NBSimulator(
    Δt::Unitful.Time,
    steps::Integer;
    t₀::Unitful.Time = 0.0 * TIME_UNIT,
    thermostat::Thermostat = NullThermostat(),
    simulator::OrdinaryDiffEqAlgorithm = VelocityVerlet(),
    potentials::Dict{Symbol,PotentialParameters} = Dict{Symbol,PotentialParameters}()
)
    NBSimulator(austrip(Δt), steps, austrip(t₀), thermostat, simulator, potentials)
end

# Extract the tuple of start_time, end_time from the simulator
time_range(simulator::NBSimulator) = (simulator.t₀, simulator.t₀ + simulator.steps * simulator.Δt)

function simulate(system::AbstractSystem, simulator::NBSimulator, potential::ArbitraryPotential)
    simulator.potentials[:custom] = InteratomicPotentialParameters(potential)
    simulate(system, simulator)
end
function simulate(system::AbstractSystem, simulator::NBSimulator, potential::LennardJonesParameters)
    simulator.potentials[:lennard_jones] = potential
    simulate(system, simulator)
end
function simulate(system::AbstractSystem, simulator::NBSimulator)
    nb_system = PotentialNBodySystem{ElementMassBody}(bodies(system), simulator.potentials)
    simulation = NBodySimulation(nb_system, time_range(simulator), nbs_boundary_conditions(system), simulator.thermostat, austrip(u"k"))
    NBSResult(run_simulation(simulation, simulator.simulator, dt = simulator.Δt))
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

# Internal struct that represents a potential in the format accepted by NBodySimulator
# Wraps the underlying InteratomicPotential with a cache of the forces for the current timestep
struct InteratomicPotentialParameters <: PotentialParameters
    potential::ArbitraryPotential
    timestep_cache::RefValue{Real}
    force_cache::RefValue{Vector{SVector{3,Real}}}
    InteratomicPotentialParameters(potential::ArbitraryPotential) = new(potential, Ref{Real}(), Ref{Vector{SVector{3,Real}}}())
end

function NBodySimulator.get_accelerating_function(parameters::InteratomicPotentialParameters, simulation::NBodySimulation)
    bounding_box = ab_bounding_box(simulation.boundary_conditions)
    boundary_conditions = ab_boundary_conditions(simulation.boundary_conditions)
    symbols = get_atomic_symbols(simulation.system)
    numbers = get_atomic_numbers(simulation.system)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin
        if !isassigned(parameters.timestep_cache) || t != parameters.timestep_cache[]
            positions = [u[:, j] for j ∈ 1:size(u, 2)] * LENGTH_UNIT
            system = FastSystem(bounding_box, boundary_conditions, positions, symbols, numbers, masses * MASS_UNIT)
            wrapped_system = DynamicSystem(system, t * TIME_UNIT)
            parameters.timestep_cache[] = t
            parameters.force_cache[] = force(wrapped_system, parameters.potential)
        end
        dv .+= parameters.force_cache[][i] / masses[i]
    end
end

# -----------------------------------------------------------------------------
# Unitful Support for Built-in LennardJonesParameters
# -----------------------------------------------------------------------------

function NBodySimulator.LennardJonesParameters(ϵ::Unitful.Energy, σ::Unitful.Length, R::Unitful.Length)
    LennardJonesParameters(austrip(ϵ), austrip(σ), austrip(R))
end