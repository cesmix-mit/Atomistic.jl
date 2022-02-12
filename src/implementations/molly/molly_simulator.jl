# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsSimulator API
# -----------------------------------------------------------------------------

@kwdef struct MollySimulator{S,T<:Unitful.Time} <: MolecularDynamicsSimulator
    Δt::T
    steps::Int
    t₀::T = zero(T)
    coupling = NoCoupling()
    parallel::Bool = false
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
    wrapper = InteratomicPotentialInteraction(potential)
    neighbor_finder = SelfNeighborFinder()
    loggers = Dict(
        "c" => CoordinateLogger(1),
        "v" => VelocityLogger(1),
        "p" => CustomPotentialLogger(1),
        "t" => TemperatureLogger(1),
        "k" => CustomKineticEnergyLogger(1)
    )
    loggers["p"].force_cache[] = InteratomicPotentials.force(system, potential) .* FORCE_UNIT  # preload initial forces
    system = System(system; general_inters = (wrapper,), neighbor_finder = neighbor_finder, loggers = loggers)
    simulate!(system, S(simulator.Δt, simulator.coupling), simulator.steps; parallel = simulator.parallel)
    MollyResult(system, simulator)
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

struct InteratomicPotentialInteraction <: GeneralInteraction
    nl_only::Bool  # must be true to use the SelfNeighborFinder hack
    potential::ArbitraryPotential
    InteratomicPotentialInteraction(potential::ArbitraryPotential) = new(true, potential)
end

@inline @inbounds function Molly.force!(fs, inter::InteratomicPotentialInteraction, s::System, i::Integer, j::Integer, force_units, weight_14::Bool = false)
    @assert i == j
    fdr = s.loggers["p"].force_cache[][i]
    Molly.check_force_units(fdr, force_units)
    fs[i] -= ustrip.(fdr)
    return nothing
end

struct SelfNeighborFinder <: AbstractNeighborFinder end

function Molly.find_neighbors(s::System, nf::SelfNeighborFinder, current_neighbors = nothing, step_n::Integer = 0; parallel::Bool = true)
    Molly.NeighborList(length(s), [(i, i, false) for i ∈ 1:length(s)])
end

struct CustomPotentialLogger{F<:Unitful.Force,E<:Unitful.Energy}
    n_steps::Int
    force_cache::Ref{Vector{SVector{3,F}}}
    energies::Vector{E}
    CustomPotentialLogger{F,E}(n_steps::Int) where {F<:Unitful.Force,E<:Unitful.Energy} = new{F,E}(n_steps, Ref{Vector{SVector{3,F}}}(), E[])
end
function CustomPotentialLogger(n_steps::Int)
    F = typeof(zeros() * FORCE_UNIT)
    E = typeof(zeros() * ENERGY_UNIT)
    CustomPotentialLogger{F,E}(n_steps)
end

function Molly.log_property!(logger::CustomPotentialLogger, s::System, neighbors = nothing, step_n::Int = 0)
    if step_n % logger.n_steps == 0
        eandf = energy_and_force(s, first(s.general_inters).potential)
        logger.force_cache[] = eandf.f .* FORCE_UNIT
        push!(logger.energies, eandf.e * ENERGY_UNIT)
    end
end

struct CustomKineticEnergyLogger{E<:Unitful.Energy}
    n_steps::Int
    energies::Vector{E}
    CustomKineticEnergyLogger{E}(n_steps::Int) where {E<:Unitful.Energy} = new{E}(n_steps, E[])
end
function CustomKineticEnergyLogger(n_steps::Int)
    E = typeof(zeros() * ENERGY_UNIT)
    CustomKineticEnergyLogger{E}(n_steps)
end

function Molly.log_property!(logger::CustomKineticEnergyLogger, s::System, neighbors = nothing, step_n::Int = 0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, Molly.kinetic_energy(s))
    end
end
