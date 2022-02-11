# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsSimulator API
# -----------------------------------------------------------------------------

@kwdef struct MollySimulator{S,T<:Unitful.Time} <: MolecularDynamicsSimulator
    Δt::T
    steps::Int
    t₀::T = zero(T)
    coupling = NoCoupling()
    parallel::Bool = false  # TODO
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
        "t" => TemperatureLogger(1),
        "k" => CustomKineticEnergyLogger(1),
        "p" => CustomPotentialLogger(1)
    )
    system = System(system; general_inters = (wrapper,), neighbor_finder = neighbor_finder, loggers = loggers)
    wrapper.system[] = system
    simulate!(system, S(simulator.Δt, simulator.coupling), simulator.steps; parallel = simulator.parallel)
    MollyResult(system, simulator)
end

# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------

struct InteratomicPotentialInteraction{F<:Unitful.Force,E<:Unitful.Energy} <: GeneralInteraction
    nl_only::Bool  # must be true to use the SelfNeighborFinder hack
    potential::ArbitraryPotential
    system::Ref{System}
    force_cache::Ref{Vector{SVector{3,F}}}
    energy_cache::Ref{E}
    function InteratomicPotentialInteraction(potential::ArbitraryPotential)
        F = typeof(zeros() * FORCE_UNIT)
        E = typeof(zeros() * ENERGY_UNIT)
        new{F,E}(true, potential, Ref{System}(), Ref{Vector{SVector{3,F}}}(), Ref{E}())
    end
end

function Molly.force(inter::InteratomicPotentialInteraction, vec_ij, coord_i, coord_j, atom_i, atom_j, box_size)
    @assert atom_i == atom_j
    if atom_i.index == 1
        eandf = energy_and_force(inter.system[], inter.potential)
        inter.force_cache[] = eandf.f .* FORCE_UNIT
        inter.energy_cache[] = eandf.e * ENERGY_UNIT
    end
    inter.force_cache[][atom_i.index]
end

struct SelfNeighborFinder <: AbstractNeighborFinder end

function Molly.find_neighbors(s::System, nf::SelfNeighborFinder, current_neighbors = nothing, step_n::Integer = 0; parallel::Bool = true)
    Molly.NeighborList(length(s), [(i, i, false) for i ∈ 1:length(s)])
end

struct CustomPotentialLogger{E<:Unitful.Energy}
    n_steps::Integer
    energies::Vector{E}
    CustomPotentialLogger{E}(n_steps::Integer) where {E<:Unitful.Energy} = new{E}(n_steps, E[])
end
function CustomPotentialLogger(n_steps::Integer)
    E = typeof(zeros() * ENERGY_UNIT)
    CustomPotentialLogger{E}(n_steps)
end

function Molly.log_property!(logger::CustomPotentialLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, first(s.general_inters).energy_cache[])
    end
end

struct CustomKineticEnergyLogger{E<:Unitful.Energy}
    n_steps::Int
    energies::Vector{E}
    CustomKineticEnergyLogger{E}(n_steps::Integer) where {E<:Unitful.Energy} = new{E}(n_steps, E[])
end
function CustomKineticEnergyLogger(n_steps::Integer)
    E = typeof(zeros() * ENERGY_UNIT)
    CustomKineticEnergyLogger{E}(n_steps)
end

function Molly.log_property!(logger::CustomKineticEnergyLogger, s::System, neighbors = nothing, step_n::Integer = 0)
    if step_n % logger.n_steps == 0
        push!(logger.energies, Molly.kinetic_energy(s))
    end
end
