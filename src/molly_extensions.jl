using LinearAlgebra
using Molly
using Parameters

export EquipartitionThermostat,
    KineticEnergyLogger,
    PotentialEnergyLogger,
    VelocityLogger

import Molly.apply_thermostat!
import Molly.log_property!
import Molly.kinetic_energy
import Molly.potential_energy
import Molly.DefaultFloat

# CUSTOM THERMOSTAT

kB = 1.380649e-23 # Boltzmann constant in J / K

@with_kw struct EquipartitionThermostat{T <: Int} <: Thermostat
    dof::T = 3
end

function apply_thermostat!(velocities, s::Simulation, thermostat::EquipartitionThermostat)
    # return [
    #     normalize(v) * √(thermostat.dof * kB * s.temperature / (s.atoms[i].mass * 1.66054e-27)) * 10^-3
    #     for (i, v) ∈ enumerate(velocities)
    # ]
    for (i, v) ∈ enumerate(velocities)
        velocities[i] = normalize(v) * √(thermostat.dof * kB * s.temperature / (s.atoms[i].mass * 1.66054e-27)) * 10^-3
    end
    return velocities
end

# CUSTOM LOGGERS

struct KineticEnergyLogger{T} <: Logger
    n_steps::Int
    energies::Vector{T}
end

KineticEnergyLogger(T::Type, n_steps::Integer) = EnergyLogger(n_steps, T[])

KineticEnergyLogger(n_steps::Integer) = KineticEnergyLogger(DefaultFloat, n_steps)

function Base.show(io::IO, el::KineticEnergyLogger)
    print(io, "KineticEnergyLogger{", eltype(el.energies), "} with n_steps ", el.n_steps, ", ", length(el.energies), " energies recorded")
end

function log_property!(logger::KineticEnergyLogger, s::Simulation, step_n::Integer)
    if step_n % logger.n_steps == 0
        push!(logger.energies, kinetic_energy(s))
    end
end

struct PotentialEnergyLogger{T} <: Logger
    n_steps::Int
    energies::Vector{T}
end

PotentialEnergyLogger(T::Type, n_steps::Integer) = PotentialEnergyLogger(n_steps, T[])

PotentialEnergyLogger(n_steps::Integer) = PotentialEnergyLogger(DefaultFloat, n_steps)

function Base.show(io::IO, el::PotentialEnergyLogger)
    print(io, "PotentialEnergyLogger{", eltype(el.energies), "} with n_steps ", el.n_steps, ", ", length(el.energies), " energies recorded")
end

function log_property!(logger::PotentialEnergyLogger, s::Simulation, step_n::Integer)
    if step_n % logger.n_steps == 0
        push!(logger.energies, potential_energy(s))
    end
end

struct VelocityLogger{T} <: Logger
    n_steps::Int
    velocities::Vector{Vector{T}}
end

function VelocityLogger(T, n_steps::Integer; dims::Integer=3)
    return VelocityLogger(n_steps, Array{SArray{Tuple{dims},T,1,dims},1}[])
end

function VelocityLogger(n_steps::Integer; dims::Integer=3)
    return VelocityLogger(DefaultFloat, n_steps, dims=dims)
end

function Base.show(io::IO, vl::VelocityLogger)
    print(io, "VelocityLogger{", eltype(eltype(vl.velocities)), "} with n_steps ", vl.n_steps, ", ", length(vl.velocities), " frames recorded for ", length(vl.velocities) > 0 ? length(first(vl.velocities)) : "?", " atoms")
end

function log_property!(logger::VelocityLogger, s::Simulation, step_n::Integer)
    if step_n % logger.n_steps == 0
        push!(logger.velocities, deepcopy(s.velocities))
    end
end