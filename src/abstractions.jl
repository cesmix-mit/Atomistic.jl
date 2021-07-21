abstract type NuclearPotentialParameters end
abstract type MolecularDynamicsParameters end

Base.@kwdef struct LJParameters <: NuclearPotentialParameters
    ϵ::Quantity
    σ::Quantity
    R::Quantity
end
