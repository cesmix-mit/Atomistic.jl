# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsSimulator API
# -----------------------------------------------------------------------------

"""
    MollySimulator{S,T<:Unitful.Time,C} <: MolecularDynamicsSimulator

A wrapper around Molly to implement the Atomistic API.

**Type parameter descriptions**
- `S` the type of simulator
- `T<:Unitful.Time` the type for the time fields
- `C` the type for the coupling field

**Field descriptions**
- `Δt::T` the time between timesteps, assumed to be in atomic units
- `steps::Integer` the number of timesteps for the simulation
- `t₀::T` the starting time of the simulation, assumed to be in atomic units;
    defaults to 0
- `coupling::C` the coupling for the simulation;
    many options are defined by `Molly`, but a user could also define a custom thermostat;
    defaults to `NoCoupling()`
- `stride::Int` the number of timesteps between logger runs
    defaults to 1
"""
struct MollySimulator{S,T<:Unitful.Time,C} <: MolecularDynamicsSimulator
    Δt::T
    steps::Int
    t₀::T
    coupling::C
    stride::Int
end
function MollySimulator{S}(Δt::T, steps::Int; t₀::Real = zero(T), coupling::C = NoCoupling(), stride::Int = 1) where {S,T<:Real,C}
    Δt, t₀ = promote(Δt * TIME_UNIT, t₀ * TIME_UNIT)
    MollySimulator{S,typeof(Δt),C}(Δt, steps, t₀, coupling, stride)
end
function MollySimulator{S}(Δt::T, steps::Int; t₀::Unitful.Time = zero(T), coupling::C = NoCoupling(), stride::Int = 1) where {S,T<:Unitful.Time,C}
    Δt, t₀ = promote(Δt, t₀)
    MollySimulator{S,T,C}(Δt, steps, t₀, coupling, stride)
end
function MollySimulator(Δt::T, steps::Int; t₀::Real = zero(T), coupling::C = NoCoupling(), stride::Int = 1) where {T<:Real,C}
    MollySimulator{Molly.VelocityVerlet}(Δt, steps; t₀ = t₀, coupling = coupling, stride = stride)
end
function MollySimulator(Δt::T, steps::Int; t₀::Unitful.Time = zero(T), coupling::C = NoCoupling(), stride::Int = 1) where {T<:Unitful.Time,C}
    MollySimulator{Molly.VelocityVerlet}(Δt, steps; t₀ = t₀, coupling = coupling, stride = stride)
end

function simulate(system::AbstractSystem{3}, simulator::MollySimulator{S}, potential::ArbitraryPotential) where {S}
    wrapper = InteratomicPotentialInter(potential)
    loggers = Dict(
        "c" => CoordinateLogger(LENGTH_TYPE, simulator.stride),
        "v" => VelocityLogger(VELOCITY_TYPE, simulator.stride),
        "p" => PotentialEnergyLogger(ENERGY_TYPE, simulator.stride),
        "k" => KineticEnergyLogger(ENERGY_TYPE, simulator.stride),
        "t" => TemperatureLogger(TEMPERATURE_TYPE, simulator.stride)
    )
    system = System(system; general_inters = (wrapper,), loggers = loggers)
    simulate!(system, S(dt = simulator.Δt, coupling = simulator.coupling), simulator.steps)
    MollyResult(system, simulator)
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

# Internal struct that represents a potential in the format accepted by Molly
# Wraps the underlying InteratomicPotential with a cache of the potential energy for the current timestep
struct InteratomicPotentialInter{P<:ArbitraryPotential,E<:Unitful.Energy}
    potential::P
    energy_cache::Ref{E}
    InteratomicPotentialInter(potential::ArbitraryPotential) = new{typeof(potential),ENERGY_TYPE}(potential, Ref{ENERGY_TYPE}())
end

function Molly.forces(inter::InteratomicPotentialInter, sys, neighbors = nothing)
    eandf = energy_and_force(sys, inter.potential)
    inter.energy_cache[] = eandf.e * sys.energy_units
    eandf.f .* sys.force_units
end

function Molly.potential_energy(inter::InteratomicPotentialInter, sys, neighbors = nothing)
    inter.energy_cache[]
end
