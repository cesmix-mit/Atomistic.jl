# Molecular Dynamics Result Abstract Interface

abstract type MolecularDynamicsResult end

function get_bodies(result::MolecularDynamicsResult, t::Integer=0)::Vector{Atom}
    throw(UnimplementedError(:get_bodies, result))
end
function get_time_range(result::MolecularDynamicsResult)::Vector{<:Real}
    throw(UnimplementedError(:get_time_range, result))
end

function plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer)::Plot
    throw(UnimplementedError(:plot_temperature!, result))
end
function plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer)::Plot
    throw(UnimplementedError(:plot_energy!, result))
end

function calculate_rdf(result::MolecularDynamicsResult, sample_fraction::Real)::Tuple{Vector{<:Real},Vector{<:Real}}
    throw(UnimplementedError(:calculate_rdf, result))
end

function plot_temperature(result::MolecularDynamicsResult, stride::Integer)
    N = length(get_bodies(result))
    p = plot(
        title="Temperature during Simulation [n = $(N)]",
        xlab="Time",
        ylab="Temperature",
    )
    plot_temperature!(p, result, stride)
end

function plot_energy(result::MolecularDynamicsResult, stride::Integer)
    N = length(get_bodies(result))
    p = plot(
        title="Energy during Simulation [n = $(N)]",
        xlab="Time",
        ylab="Energy",
        legend=:right
    ) 
    plot_energy!(p, result, stride)
end

function plot_rdf(result::MolecularDynamicsResult, σ::Real, sample_fraction::Real=0.1)
    plot_rdf(result, σ * u"bohr", sample_fraction)
end
function plot_rdf(result::MolecularDynamicsResult, σ::Quantity, sample_fraction::Real=0.1)
    N = length(get_bodies(result))
    T = length(get_time_range(result)) - 1
    rs, grf = calculate_rdf(result, sample_fraction)
    plot(
        title="Radial Distribution Function [n = $(N)] [T = $(T)]",
        xlab="Distance r/σ",
        ylab="Radial Distribution g(r)",
        legend=false
    )
    plot!(
        rs * u"bohr" / auconvert(σ),
        grf
    )
end