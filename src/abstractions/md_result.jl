# Molecular Dynamics Result Abstract Interface

abstract type MolecularDynamicsResult end

function get_system(result::MolecularDynamicsResult, t::Integer=0)::AbstractSystem
    throw(UnimplementedError(:get_system, result))
end
function get_time_range(result::MolecularDynamicsResult)::Vector{<:Real}
    throw(UnimplementedError(:get_time_range, result))
end

function temperature(result::MolecularDynamicsResult, time::Real)::Real
    throw(UnimplementedError(:temperature, result))
end
reference_temperature(result::MolecularDynamicsResult)::Union{Real,Missing} = missing

function kinetic_energy(result::MolecularDynamicsResult, time::Real)::Real
    throw(UnimplementedError(:kinetic_energy, result))
end
function potential_energy(result::MolecularDynamicsResult, time::Real)::Real
    throw(UnimplementedError(:potential_energy, result))
end
function total_energy(result::MolecularDynamicsResult, time::Real)::Real
    throw(UnimplementedError(:total_energy, result))
end

function rdf(result::MolecularDynamicsResult, sample_fraction::Float64=1.0)::Tuple{Vector{<:Real},Vector{<:Real}}
    throw(UnimplementedError(:rdf, result))
end