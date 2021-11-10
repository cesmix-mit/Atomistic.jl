# functions for plotting results

function plot_temperature(result::MolecularDynamicsResult, stride::Integer)
    N = length(get_system(result))
    p = plot(
        title="Temperature during Simulation [n = $(N)]",
        xlab="Time",
        ylab="Temperature",
    )
    plot_temperature!(p, result, stride, true)
end

function plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer, first_plot::Bool=false)
    time_range = [(t, auconvert(u"ps", time)) for (t, time) ∈ enumerate(get_time_range(result)) if (t - 1) % stride == 0]
    if (!first_plot)
        vline!(
            p,
            [time_range[1][2]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [auconvert(u"K", temperature(result, austrip(t))) for (t, time) ∈ time_range],
        label=first_plot ? "Simulation Temperature" : nothing,
        color=1,
    )
    reference_temp = reference_temperature(result)
    if (!ismissing(reference_temp))
        plot!(
            p,
            [time for (t, time) ∈ time_range],
            [auconvert(u"K", reference_temp) for (t, time) ∈ time_range],
            label=first_plot ? "Reference Temperature" : nothing,
            color=2,
            linestyle=:dash,
            lw=2
        )
    end
    p
end

function plot_energy(result::MolecularDynamicsResult, stride::Integer)
    N = length(get_system(result))
    p = plot(
        title="Energy during Simulation [n = $(N)]",
        xlab="Time",
        ylab="Energy",
        legend=:right
    ) 
    plot_energy!(p, result, stride, true)
end

function plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer, first_plot::Bool=false)
    time_range = [(t, auconvert(u"ps", time)) for (t, time) ∈ enumerate(get_time_range(result)) if (t - 1) % stride == 0]
    if (!first_plot)
        vline!(
            p,
            [time_range[1][2]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [kinetic_energy(result, austrip(t))u"hartree" for (t, time) ∈ time_range],
        label=first_plot ? "Kinetic Energy" : nothing,
        color=2
    )
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [potential_energy(result, austrip(t))u"hartree" for (t, time) ∈ time_range],
        label=first_plot ? "Potential Energy" : nothing,
        color=1
    )
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [total_energy(result, austrip(t))u"hartree" for (t, time) ∈ time_range],
        label=first_plot ? "Total Energy" : nothing,
        color=3
    )
end

function plot_rdf(result::MolecularDynamicsResult, σ::Real, sample_fraction::Float64=1.0)
    plot_rdf(result, σ * u"bohr", sample_fraction)
end
function plot_rdf(result::MolecularDynamicsResult, σ::Unitful.Length, sample_fraction::Float64=1.0)
    @assert 0 < sample_fraction ≤ 1
    N = length(get_system(result))
    T = length(get_time_range(result)) - 1
    rs, grf = rdf(result, sample_fraction)
    @assert length(rs) == length(grf)
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
