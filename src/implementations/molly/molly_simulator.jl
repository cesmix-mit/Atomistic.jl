# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsSimulator API
# -----------------------------------------------------------------------------

@kwdef struct MollySimulator{S,T<:Unitful.Time} <: MolecularDynamicsSimulator
    Δt::T
    steps::Int
    t₀::T = zero(T)
    coupling = NoCoupling()
end
function MollySimulator{S}(Δt::T, steps::Int; t₀::Real = zero(T), kwargs...) where {S,T<:Real}
    Δt, t₀ = promote(Δt * TIME_UNIT, t₀ * TIME_UNIT)
    MollySimulator{S,typeof(Δt)}(; Δt = Δt, steps = steps, t₀ = t₀, kwargs...)
end
function MollySimulator{S}(Δt::T, steps::Integer; t₀::Unitful.Time = zero(T), kwargs...) where {S,T<:Unitful.Time}
    Δt, t₀ = promote(Δt, t₀)
    MollySimulator{S,T}(; Δt = Δt, steps = steps, t₀ = t₀, kwargs...)
end
function MollySimulator(Δt::T, steps::Int; t₀::Real = zero(T), kwargs...) where {T<:Real}
    Δt, t₀ = promote(Δt * TIME_UNIT, t₀ * TIME_UNIT)
    MollySimulator{Molly.VelocityVerlet,typeof(Δt)}(; Δt = Δt, steps = steps, t₀ = t₀, kwargs...)
end
function MollySimulator(Δt::T, steps::Integer; t₀::Unitful.Time = zero(T), kwargs...) where {T<:Unitful.Time}
    Δt, t₀ = promote(Δt, t₀)
    MollySimulator{Molly.VelocityVerlet,T}(; Δt = Δt, steps = steps, t₀ = t₀, kwargs...)
end

function simulate(system::AbstractSystem{3}, simulator::MollySimulator{S}, potential::ArbitraryPotential) where {S}
    wrapper = InteratomicPotentialInter(potential)
    loggers = Dict(
        "c" => CoordinateLogger(typeof(zeros() * LENGTH_UNIT), 1),
        "v" => VelocityLogger(typeof(zeros() * VELOCITY_UNIT), 1),
        "p" => PotentialEnergyLogger(typeof(zeros() * ENERGY_UNIT), 1),
        "k" => KineticEnergyLogger(typeof(zeros() * ENERGY_UNIT), 1),
        "t" => TemperatureLogger(typeof(zeros() * TEMPERATURE_UNIT), 1)
    )
    system = System(system; general_inters = (wrapper,), loggers = loggers)
    simulate!(system, S(simulator.Δt, simulator.coupling), simulator.steps)
    MollyResult(system, simulator)
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

struct InteratomicPotentialInter{E<:Unitful.Energy}
    potential::ArbitraryPotential
    energy_cache::Ref{E}
    function InteratomicPotentialInter(potential::ArbitraryPotential)
        E = typeof(zeros() * ENERGY_UNIT)
        new{E}(potential, Ref{E}())
    end
end

function Molly.forces(inter::InteratomicPotentialInter, sys, neighbors = nothing)
    eandf = energy_and_force(sys, inter.potential)
    inter.energy_cache[] = eandf.e * sys.energy_units
    eandf.f .* sys.force_units
end

function Molly.potential_energy(inter::InteratomicPotentialInter, sys, neighbors = nothing)
    inter.energy_cache[]
end
