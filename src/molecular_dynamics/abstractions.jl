# Molecular Dynamics Abstractions

abstract type MolecularDynamicsParameters end

abstract type MolecularDynamicsResult end

function simulate(bodies::AbstractVector{<:MassBody}, parameters::MolecularDynamicsResult, nuclear_potential_parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:simulate, parameters))
end

function get_bodies(result::MolecularDynamicsResult, t::Integer=0)
    throw(UnimplementedError(:final_bodies, result))
end

function get_time_range(result::MolecularDynamicsResult)
    throw(UnimplementedError(:get_time_range, result))
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

function plot_temperature!(p::Plots.Plot, result::MolecularDynamicsResult, stride::Integer)
    throw(UnimplementedError(:plot_temperature!, result))
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

function plot_energy!(p::Plots.Plot, result::MolecularDynamicsResult, stride::Integer)
    throw(UnimplementedError(:plot_energy!, result))
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

function calculate_rdf(result::MolecularDynamicsResult, sample_fraction::Real)
    throw(UnimplementedError(:calculate_rdf, result))
end
