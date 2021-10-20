# Molecular Dynamics Result Abstract Interface

abstract type MolecularDynamicsResult end

function get_bodies(result::MolecularDynamicsResult, t::Integer=0)::Vector{Atom}
    throw(UnimplementedError(:get_bodies, result))
end
function get_time_range(result::MolecularDynamicsResult)::Vector{<:Real}
    throw(UnimplementedError(:get_time_range, result))
end

function temperature(result::MolecularDynamicsResult, time::Real)::Real
    throw(UnimplementedError(:temperature, result))
end
function reference_temperature(result::MolecularDynamicsResult)::Union{Real,Nothing}
    nothing
end

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