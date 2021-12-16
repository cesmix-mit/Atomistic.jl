# Integrations with NBodySimulator.jl

struct ElementMassBody{cType<:Real,mType<:Real} <: Body
    r::SVector{3,cType}
    v::SVector{3,cType}
    m::mType
    s::Symbol
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, e::Element)
    ElementMassBody{Float64,Float64}(austrip.(r), austrip.(v), austrip(e.atomic_mass), Symbol(e.symbol))
end

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

function get_symbols(system::PotentialNBodySystem{ElementMassBody})
    [system.bodies[i].s for i ∈ 1:length(system.bodies)]
end

function bodies(system::AbstractSystem)
    [ElementMassBody(position(a), velocity(a), species(a)) for a ∈ system]
end

function nbody_boundary_conditions(system::AbstractSystem)
    # TODO: support more boundary conditions
    box = bounding_box(system)
    @assert box[1][1] == box[2][2] == box[3][3]
    @assert boundary_conditions(system) == [Periodic(), Periodic(), Periodic()]
    CubicPeriodicBoundaryConditions(austrip(box[1][1]))
end

# TODO: support more boundary conditions
function DynamicAtom(b::ElementMassBody, boundary_conditions::CubicPeriodicBoundaryConditions)
    DynamicAtom(b, boundary_conditions.L * u"bohr")
end
function DynamicAtom(b::ElementMassBody, box_size::Unitful.Length)
    DynamicAtom(mod.(b.r .* u"bohr", box_size), b.v .* u"bohr * hartree / ħ_au", elements[b.s])
end

# TODO: support more boundary conditions
function DynamicSystem(bodies::Vector{<:ElementMassBody}, boundary_conditions::CubicPeriodicBoundaryConditions, time::Real = 0.0)
    DynamicSystem(bodies, boundary_conditions.L * u"bohr", time * u"ħ_au / hartree")
end
function DynamicSystem(bodies::Vector{<:ElementMassBody}, box_size::Unitful.Length, time::Unitful.Time = 0.0u"s")
    DynamicSystem(DynamicAtom.(bodies, box_size), box_size, time)
end

@kwdef struct InteratomicPotentialParameters <: PotentialParameters
    potential::ArbitraryPotential
    timestep_cache::RefValue{Real} = Ref{Real}()
    force_cache::RefValue{Vector{SVector{3,Real}}} = Ref{Vector{SVector{3,Real}}}()
end

function NBodySimulator.get_accelerating_function(parameters::InteratomicPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    symbols = get_symbols(simulation.system)
    (dv, u, v, t, i) -> begin
        if !isassigned(parameters.timestep_cache) || t != parameters.timestep_cache[]
            particles = [ElementMassBody(SVector{3}(u[:, j]), SVector{3}(v[:, j]), masses[j], symbols[j]) for j ∈ 1:length(symbols)]
            system = DynamicSystem(particles, simulation.boundary_conditions, t)
            parameters.timestep_cache[] = t
            parameters.force_cache[] = force(system, parameters.potential)
        end
        dv .+= parameters.force_cache[][i] / masses[i]
    end
end

@kwdef struct LJPotential
    ϵ::Real
    σ::Real
    R::Real
end
function LJPotential(ϵ::Unitful.Energy, σ::Unitful.Length, R::Unitful.Length)
    LJPotential(austrip(ϵ), austrip(σ), austrip(R))
end

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
    Δt::Real
    steps::Integer
    t₀::Real = 0.0
    thermostat::Thermostat = NullThermostat()
    simulator::OrdinaryDiffEqAlgorithm = VelocityVerlet()
    potentials::Dict{Symbol,PotentialParameters} = Dict{Symbol,PotentialParameters}()
end
function NBSimulator(Δt::Unitful.Time,
    steps::Integer,
    t₀::Unitful.Time = 0.0u"s",
    thermostat::Thermostat = NullThermostat(),
    simulator::OrdinaryDiffEqAlgorithm = VelocityVerlet(),
    potentials::Dict{Symbol,PotentialParameters} = Dict{Symbol,PotentialParameters}())
    NBSimulator(austrip(Δt), steps, austrip(t₀), thermostat, simulator, potentials)
end

function simulate(system::AbstractSystem, simulator::NBSimulator, potential::ArbitraryPotential)
    simulator.potentials[:custom] = InteratomicPotentialParameters(potential = potential)
    simulate(system, simulator)
end
function simulate(system::AbstractSystem, simulator::NBSimulator, potential::LJPotential)
    simulator.potentials[:lennard_jones] = LennardJonesParameters(potential.ϵ, potential.σ, potential.R)
    simulate(system, simulator)
end
function simulate(system::AbstractSystem, simulator::NBSimulator)
    nb_system = PotentialNBodySystem{ElementMassBody}(bodies(system), simulator.potentials)
    simulation = NBodySimulation(nb_system, (simulator.t₀, simulator.t₀ + simulator.steps * simulator.Δt), nbody_boundary_conditions(system), simulator.thermostat, 1.0)
    NBSResult(run_simulation(simulation, simulator.simulator, dt = simulator.Δt))
end

"""
    NBSResult <: MolecularDynamicsResult

The result generating from running a `MolecularDynamicsSimulator`.

**Field descriptions**
- `result::SimulationResult` the standard simulation result from `NBodySimulator`
"""
struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
end

function get_system(result::NBSResult, t::Integer = 0)
    sr = result.result
    time = get_time_range(result)[t > 0 ? t : end]
    positions = get_position(sr, time)
    velocities = get_velocity(sr, time)
    masses = get_masses(sr.simulation.system)
    symbols = get_symbols(sr.simulation.system)
    particles = [ElementMassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i], symbols[i]) for i ∈ 1:length(symbols)]
    DynamicSystem(particles, sr.simulation.boundary_conditions, time)
end

function get_time_range(result::NBSResult)
    result.result.solution.t
end

function temperature(result::NBSResult, t::Integer = 0)
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.temperature(result.result, time)
end

function reference_temperature(result::NBSResult)
    result.result.simulation.thermostat isa NullThermostat ? missing : result.result.simulation.thermostat.T
end

function kinetic_energy(result::NBSResult, t::Integer = 0)
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.kinetic_energy(result.result, time)
end

function potential_energy(result::NBSResult, t::Integer = 0)
    potentials = result.result.simulation.system.potentials
    # https://github.com/SciML/NBodySimulator.jl/issues/44
    if :custom ∈ keys(potentials)
        return InteratomicPotentials.potential_energy(get_system(result, t), potentials[:custom].potential)
    end
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.potential_energy(result.result, time)
end

function rdf(result::NBSResult, sample_fraction::Float64 = 1.0)
    @assert 0 < sample_fraction ≤ 1
    sr = result.result
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    trange = get_time_range(result)[end-floor(Int, length(sr.solution.t) * sample_fraction)+1:end]

    maxbin = 1000
    dr = pbc.L / 2 / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for i ∈ 1:n
            ri = SVector{3}(cc[:, i])
            for j ∈ i+1:n
                rj = SVector{3}(cc[:, j])
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
