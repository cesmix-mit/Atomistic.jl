# This file will be superceded by Pontentials.jl

# Nuclear Potential Abstractions

abstract type NuclearPotentialParameters end

function forces(state::AtomicConfiguration, parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:forces, parameters))
end

function potential_energy(state::AtomicConfiguration, parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:potential_energy, parameters))
end

Base.@kwdef struct LJParameters <: NuclearPotentialParameters
    ϵ::Quantity
    σ::Quantity
    R::Quantity
end
