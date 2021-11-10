# Integrations with NBodySimulator.jl

struct ElementMassBody{cType <: Real,mType <: Real} <: Body
    NBodySimulator.@position_velocity_mass
    e::ChemicalElement
end
function ElementMassBody(r::SVector{3,<:Unitful.Length}, v::SVector{3,<:Unitful.Velocity}, e::ChemicalElement)
    ElementMassBody(austrip.(r), austrip.(v), austrip(atomic_mass(e)), e)
end

function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, e::ChemicalElement, L::Unitful.Length, reference_temp::Unitful.Temperature; rng=MersenneTwister(n))
    average_velocity = √(u"k" * reference_temp / atomic_mass(e))
    generate_bodies_in_cell_nodes(n, e, average_velocity, L)
end
function NBodySimulator.generate_bodies_in_cell_nodes(n::Integer, e::ChemicalElement, average_velocity::Unitful.Velocity, L::Unitful.Length; rng=MersenneTwister(n))
    velocities = average_velocity * randn(rng, Float64, (3, n))
    bodies = ElementMassBody[]

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x ∈ dL / 2:dL:L, y ∈ dL / 2:dL:L, z ∈ dL / 2:dL:L
        if count > n break end
        push!(bodies, ElementMassBody(SVector(x, y, z), SVector{3}(velocities[:, count]), e))
        count += 1
    end
    return bodies
end

function get_elements(system::PotentialNBodySystem{ElementMassBody})
    [system.bodies[i].e for i ∈ 1:length(system.bodies)]
end

function bodies(system::AbstractSystem)
    [ElementMassBody(AtomsBase.position(a), velocity(a), element(a)) for a ∈ system]
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
    DynamicAtom(mod.(b.r .* u"bohr", box_size), b.v .* u"bohr * hartree / ħ_au", b.e)
end

# TODO: support more boundary conditions
function AtomsBase.SimpleSystem(bodies::Vector{<:ElementMassBody}, boundary_conditions::CubicPeriodicBoundaryConditions)
    SimpleSystem(bodies, boundary_conditions.L * u"bohr")
end
function AtomsBase.SimpleSystem(bodies::Vector{<:ElementMassBody}, box_size::Unitful.Length)
    SimpleSystem(DynamicAtom.(bodies, box_size), box_size)
end

@kwdef struct InteratomicPotentialParameters <: PotentialParameters
    potential::ArbitraryPotential
    timestep_cache::RefValue{Real} = Ref{Real}()
    force_cache::RefValue{Vector{SVector{3,Real}}} = Ref{Vector{SVector{3,Real}}}()
end

function NBodySimulator.get_accelerating_function(parameters::InteratomicPotentialParameters, simulation::NBodySimulation)
    elements = get_elements(simulation.system)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin
        if !isassigned(parameters.timestep_cache) || t != parameters.timestep_cache[]
            particles = [ElementMassBody(SVector{3}(u[:, j]), SVector{3}(v[:, j]), masses[j], elements[j]) for j ∈ 1:length(elements)]
            system = SimpleSystem(particles, simulation.boundary_conditions)
            parameters.timestep_cache[] = t
            parameters.force_cache[] = force(system, parameters.potential)
        end
        dv .+= parameters.force_cache[][i] / masses[i]
    end
end

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
                     t₀::Unitful.Time=0.0u"s",
                     thermostat::Thermostat=NullThermostat(),
                     simulator::OrdinaryDiffEqAlgorithm=VelocityVerlet(),
                     potentials::Dict{Symbol,PotentialParameters}=Dict{Symbol,PotentialParameters}())
    NBSimulator(austrip(Δt), steps, austrip(t₀), thermostat, simulator, potentials)
end
    
function simulate(system::AbstractSystem, simulator::NBSimulator, potential::ArbitraryPotential)
    simulator.potentials[:custom] = InteratomicPotentialParameters(potential=potential)
    simulate(system, simulator)
end

function simulate(system::AbstractSystem, simulator::NBSimulator)
    nb_system = PotentialNBodySystem{ElementMassBody}(bodies(system), simulator.potentials)
    simulation = NBodySimulation(nb_system, (simulator.t₀, simulator.t₀ + simulator.steps * simulator.Δt), nbody_boundary_conditions(system), simulator.thermostat, 1.0)
    NBSResult(run_simulation(simulation, simulator.simulator, dt=simulator.Δt))
end

@kwdef struct LJPotential
    ϵ::Real
    σ::Real
    R::Real
end
function LJPotential(ϵ::Unitful.Energy, σ::Unitful.Length, R::Unitful.Length)
    LJPotential(austrip(ϵ), austrip(σ), austrip(R))
end

function simulate(system::AbstractSystem, simulator::NBSimulator, potential::LJPotential)
    simulator.potentials[:lennard_jones] = LennardJonesParameters(potential.ϵ, potential.σ, potential.R)
    simulate(system, simulator)
end

struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
end

function get_system(result::NBSResult, t::Integer=0)
    sr = result.result
    time = get_time_range(result)[t > 0 ? t : end]
    positions = get_position(sr, time)
    velocities = get_velocity(sr, time)
    masses = get_masses(sr.simulation.system)
    elements = get_elements(sr.simulation.system)
    particles = [ElementMassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i], elements[i]) for i ∈ 1:length(elements)]
    SimpleSystem(particles, sr.simulation.boundary_conditions)
end

function get_time_range(result::NBSResult)
    result.result.solution.t
end

function temperature(result::NBSResult, t::Integer=0)
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.temperature(result.result, time)
end

function reference_temperature(result::NBSResult)
    result.result.simulation.thermostat isa NullThermostat ? missing : result.result.simulation.thermostat.T
end

function kinetic_energy(result::NBSResult, t::Integer=0)
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.kinetic_energy(result.result, time)
end

function potential_energy(result::NBSResult, t::Integer=0)
    potentials = result.result.simulation.system.potentials
    # https://github.com/SciML/NBodySimulator.jl/issues/44
    if :custom ∈ keys(potentials)
        return InteratomicPotentials.potential_energy(get_system(result, t), potentials[:custom].potential)
    end
    time = get_time_range(result)[t > 0 ? t : end]
    NBodySimulator.potential_energy(result.result, time)
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
            ri = SVector{3}(cc[:, i])
            for j ∈ i + 1:n
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
