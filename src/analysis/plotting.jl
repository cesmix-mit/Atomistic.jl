# functions for plotting results

function plot_temperature(result::MolecularDynamicsResult, stride::Integer)
    N = length(get_system(result))
    p = plot(
        title="Temperature during Simulation [n = $(N)]",
        xlab="Time",
        ylab="Temperature",
    )
    plot_temperature!(p, result, stride)
end

function plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer)
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(get_time_range(result)) if (i - 1) % stride == 0]
    firstPlot = austrip(time_range[1]) == 0
    if (!firstPlot)
        vline!(
            p,
            [time_range[1]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        time_range,
        t -> auconvert(u"K", temperature(result, austrip(t))),
        label=firstPlot ? "Simulation Temperature" : nothing,
        color=1,
    )
    reference_temp = reference_temperature(result)
    if (!ismissing(reference_temp))
        plot!(
            p,
            time_range,
            t -> auconvert(u"K", reference_temp),
            label=firstPlot ? "Reference Temperature" : nothing,
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
    plot_energy!(p, result, stride)
end

function plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer)
    time_range = [auconvert(u"ps", t) for (i, t) ∈ enumerate(get_time_range(result)) if (i - 1) % stride == 0]
    firstPlot = austrip(time_range[1]) == 0
    if (!firstPlot)
        vline!(
            p,
            [time_range[1]],
            label=false,
            color=:black,
            linestyle=:dot,
            lw=2
        )
    end
    plot!(
        p,
        time_range,
        t -> kinetic_energy(result, austrip(t))u"hartree",
        label=firstPlot ? "Kinetic Energy" : nothing,
        color=2
    )
    plot!(
        p,
        time_range,
        t -> potential_energy(result, austrip(t))u"hartree",
        label=firstPlot ? "Potential Energy" : nothing,
        color=1
    )
    plot!(
        p,
        time_range,
        t -> total_energy(result, austrip(t))u"hartree",
label=firstPlot ? "Total Energy" : nothing,
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
