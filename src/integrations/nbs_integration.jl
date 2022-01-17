# Integrations with NBodySimulator.jl

# -----------------------------------------------------------------------------
# Integration with AtomsBase 
# -----------------------------------------------------------------------------

# Internal type used to represent a MassBody and retain the element metadata
struct ElementMassBody{cType<:Real,mType<:Real} <: Body
    r::SVector{3,cType} # in LENGTH_UNIT
    v::SVector{3,cType} # in VELOCITY_UNIT
    m::mType            # in MASS_UNIT
    s::Symbol
    n::Int
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, e::Element)
    ElementMassBody{Float64,Float64}(austrip.(r), austrip.(v), austrip(e.atomic_mass), Symbol(e.symbol), e.number)
end

# Convert AtomsBase Atom to NBodySimulator body
function ElementMassBody(atom::Atom)
    ElementMassBody{Float64,Float64}(austrip.(position(atom)), austrip.(velocity(atom)), austrip(atomic_mass(atom)), AtomsBase.atomic_symbol(atom), atomic_number(atom))
end
# Convert NBodySimulator body to AtomsBase Atom
# TODO: support more boundary conditions
function AtomsBase.Atom(b::ElementMassBody, boundary_conditions::CubicPeriodicBoundaryConditions)
    Atom(b.s, mod.(b.r, boundary_conditions.L) .* LENGTH_UNIT, b.v .* VELOCITY_UNIT)
end

# Convert AtomsBase boundary conditions to NBodySimulator boundary conditions
function nbs_boundary_conditions(system::AbstractSystem)
    # TODO: support more boundary conditions
    box = bounding_box(system)
    @assert box[1][1] == box[2][2] == box[3][3]
    @assert boundary_conditions(system) == [Periodic(), Periodic(), Periodic()]
    CubicPeriodicBoundaryConditions(austrip(box[1][1]))
end
# Convert NBodySimulator boundary conditions to AtomsBase boundary conditions
# TODO: support more boundary conditions
ab_boundary_conditions(::CubicPeriodicBoundaryConditions) = [Periodic(), Periodic(), Periodic()]

# Convert NBodySimulator boundary conditions to AtomsBase bounding box
# TODO: support more boundary conditions
function ab_bounding_box(boundary_conditions::CubicPeriodicBoundaryConditions)
    box_size = boundary_conditions.L * LENGTH_UNIT
    z = zero(typeof(box_size))
    [[box_size, z, z], [z, box_size, z], [z, z, box_size]]
end

# Convert AtomsBase AbstractSystem to Vector of NBodySimulator bodies
function bodies(system::AbstractSystem)
    ElementMassBody.(system)
end
# Convert Vector of NBodySimulator bodies to AtomsBase FlexibleSystem
function AtomsBase.FlexibleSystem(bodies::AbstractVector{<:ElementMassBody}, boundary_conditions::BoundaryConditions)
    particles = Fix2(Atom, boundary_conditions).(bodies)
    FlexibleSystem(particles, ab_bounding_box(boundary_conditions), ab_boundary_conditions(boundary_conditions))
end
# Convert Vector of NBodySimulator bodies to AtomsBase FastSystem
function AtomsBase.FastSystem(bodies::AbstractVector{<:ElementMassBody}, boundary_conditions::BoundaryConditions)
    particles = Fix2(Atom, boundary_conditions).(bodies)
    FastSystem(particles, ab_bounding_box(boundary_conditions), ab_boundary_conditions(boundary_conditions))
end

# Extract the atomic symbol of each body in the system
get_atomic_symbols(system::PotentialNBodySystem{ElementMassBody}) = [b.s for b ∈ system.bodies]
# Extract the atomic number of each body in the system
get_atomic_numbers(system::PotentialNBodySystem{ElementMassBody}) = [b.n for b ∈ system.bodies]

# -----------------------------------------------------------------------------
# Convenience methods for generating starting configurations with element data
# -----------------------------------------------------------------------------

function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, symbol::Union{Integer,AbstractString,Symbol}, L::Unitful.Length, reference_temp::Unitful.Temperature; rng = MersenneTwister(n))
    e = elements[symbol]
    average_velocity = √(u"k" * reference_temp / e.atomic_mass)
    generate_bodies_in_cell_nodes(n, symbol, average_velocity, L, rng = rng)
end
function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, symbol::Union{Integer,AbstractString,Symbol}, average_velocity::Unitful.Velocity, L::Unitful.Length; rng = MersenneTwister(n))
    velocities = average_velocity * randn(rng, Float64, (3, n))
    e = elements[symbol]
    bodies = ElementMassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x ∈ dL/2:dL:L, y ∈ dL/2:dL:L, z ∈ dL/2:dL:L
        if count > n
            break
        end
        push!(bodies, ElementMassBody(SVector(x, y, z), SVector{3}(velocities[:, count]), e))
        count += 1
    end
    return bodies
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

# Struct to specify use of internal NBodySimulator version of Lennard Jones instead of the InteratomicPotentials version
struct LJPotential
    ϵ::Real # in ENERGY_UNIT
    σ::Real # in LENGTH_UNIT
    R::Real # in LENGTH_UNIT
end
LJPotential(ϵ::Unitful.Energy, σ::Unitful.Length, R::Unitful.Length) = LJPotential(austrip(ϵ), austrip(σ), austrip(R))

# -----------------------------------------------------------------------------
# Implementation of Atomistic API
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
function simulate(system::AbstractSystem, simulator::NBSimulator, potential::LJPotential)
    simulator.potentials[:lennard_jones] = LennardJonesParameters(potential.ϵ, potential.σ, potential.R)
    simulate(system, simulator)
end
function simulate(system::AbstractSystem, simulator::NBSimulator)
    nb_system = PotentialNBodySystem{ElementMassBody}(bodies(system), simulator.potentials)
    simulation = NBodySimulation(nb_system, time_range(simulator), nbs_boundary_conditions(system), simulator.thermostat, austrip(u"k"))
    NBSResult(run_simulation(simulation, simulator.simulator, dt = simulator.Δt))
end

"""
    NBSResult <: MolecularDynamicsResult

The result generated from running an `NBSimulator`.

**Field descriptions**
- `result::SimulationResult` the standard simulation result from `NBodySimulator`
"""
struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
end

function get_system(result::NBSResult, t::Integer = 0)
    sr = result.result
    simulation = sr.simulation
    system = simulation.system
    time = sr.solution.t[t > 0 ? t : end]
    positions = get_position(sr, time)
    velocities = get_velocity(sr, time)
    masses = get_masses(system)
    symbols = get_atomic_symbols(system)
    numbers = get_atomic_numbers(system)
    particles = [ElementMassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i], symbols[i], numbers[i]) for i ∈ 1:length(system.bodies)]
    DynamicSystem(FlexibleSystem(particles, simulation.boundary_conditions), time * TIME_UNIT)
end

function get_time_range(result::NBSResult)
    result.result.solution.t * TIME_UNIT
end

function temperature(result::NBSResult, t::Integer = 0)
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.temperature(result.result, time) * TEMPERATURE_UNIT
end

function reference_temperature(result::NBSResult)
    thermostat = result.result.simulation.thermostat
    thermostat isa NullThermostat ? missing : thermostat.T * TEMPERATURE_UNIT
end

function kinetic_energy(result::NBSResult, t::Integer = 0)
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.kinetic_energy(result.result, time) * ENERGY_UNIT
end

function potential_energy(result::NBSResult, t::Integer = 0)
    potentials = result.result.simulation.system.potentials
    # https://github.com/SciML/NBodySimulator.jl/issues/44
    if :custom ∈ keys(potentials)
        return InteratomicPotentials.potential_energy(get_system(result, t), potentials[:custom].potential) * ENERGY_UNIT
    end
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.potential_energy(result.result, time) * ENERGY_UNIT
end

function rdf(result::NBSResult, sample_fraction::Float64 = 1.0)
    @assert 0 < sample_fraction ≤ 1
    sr = result.result
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    trange = sr.solution.t[end-floor(Int, length(sr.solution.t) * sample_fraction)+1:end]

    maxbin = 1000
    dr = pbc.L / 2 / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for i ∈ 1:n
            ri = SVector{3}(cc[:, i])
            for j ∈ i+1:n
                rj = SVector{3}(cc[:, j])
                # TODO: potential for major performance improvements here
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
