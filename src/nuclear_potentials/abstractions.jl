# Nuclear Potential Abstractions

abstract type NuclearPotentialParameters end

function generate_forces(bodies::AbstractVector{<:MassBody}, parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:generate_forces, parameters))
end

Base.@kwdef struct LJParameters <: NuclearPotentialParameters
    ϵ::Quantity
    σ::Quantity
    R::Quantity
end
