# Functions for plotting molecular dynamics simulation results

PLOTTING_TIME_UNIT = u"ps"
PLOTTING_TEMPERATURE_UNIT = u"K"
PLOTTING_ENERGY_UNIT = u"hartree"

# Extract (timestep, time) tuples at stride intervals from a simulation result
function plotting_time_range(result::MolecularDynamicsResult, stride::Integer; time_unit = PLOTTING_TIME_UNIT)
    [(t, uconvert(time_unit, time)) for (t, time) ∈ enumerate(get_time_range(result)) if (t - 1) % stride == 0]
end

"""
    plot_temperature(result::MolecularDynamicsResult, stride::Integer; time_unit = u"ps", temperature_unit = u"K")::Plot

Plot the temperature of a `MolecularDynamicsResult` against time, sampling every `stride` points.
Converts time and temperature quantities to the specified units.
"""
function plot_temperature(result::MolecularDynamicsResult, stride::Integer; time_unit = PLOTTING_TIME_UNIT, temperature_unit = PLOTTING_TEMPERATURE_UNIT)
    N = get_num_bodies(result)
    p = plot(
        title = "Temperature during Simulation [n = $(N)]",
        xlab = "Time",
        ylab = "Temperature",
    )
    plot_temperature!(p, result, stride; first_plot = true, time_unit = time_unit, temperature_unit = temperature_unit)
end
"""
    plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer; first_plot::Bool = false; first_plot::Bool = false, time_unit = u"ps", temperature_unit = u"K")::Plot

Plot the temperature of a `MolecularDynamicsResult` against time, sampling every `stride` points.
Converts time and temperature quantities to the specified units.
Adds the new line to an existing plot and update the legend only if it is the first plot.
If it is not the first plot, adds a vertical line to differentiate the segments of the simulation in the plot.
"""
function plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer; first_plot::Bool = false, time_unit = PLOTTING_TIME_UNIT, temperature_unit = PLOTTING_TEMPERATURE_UNIT)
    time_range = plotting_time_range(result, stride; time_unit = time_unit)
    if (!first_plot)
        vline!(
            p,
            [time_range[1][2]],
            label = false,
            color = :black,
            linestyle = :dot,
            lw = 2
        )
    end
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [uconvert(temperature_unit, temperature(result, t)) for (t, time) ∈ time_range],
        label = first_plot ? "Simulation Temperature" : nothing,
        color = 1,
    )
    reference_temp = uconvert(temperature_unit, reference_temperature(result))
    if (!ismissing(reference_temp))
        plot!(
            p,
            [time for (t, time) ∈ time_range],
            [reference_temp for (t, time) ∈ time_range],
            label = first_plot ? "Reference Temperature" : nothing,
            color = 2,
            linestyle = :dash,
            lw = 2
        )
    end
    p
end

"""
    plot_energy(result::MolecularDynamicsResult, stride::Integer; time_unit = u"ps", energy_unit = u"hartree")::Plot

Plot the kinetic, potential, and total energy of a `MolecularDynamicsResult` against time, sampling every `stride` points.
Converts time and energy quantities to the specified units.
"""
function plot_energy(result::MolecularDynamicsResult, stride::Integer; time_unit = PLOTTING_TIME_UNIT, energy_unit = PLOTTING_ENERGY_UNIT)
    N = get_num_bodies(result)
    p = plot(
        title = "Energy during Simulation [n = $(N)]",
        xlab = "Time",
        ylab = "Energy",
        legend = :right
    )
    plot_energy!(p, result, stride; first_plot = true, time_unit = time_unit, energy_unit = energy_unit)
end
"""
    plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer; first_plot::Bool = false, time_unit = u"ps", energy_unit = u"hartree")::Plot

Plot the kinetic, potential, and total energy of a `MolecularDynamicsResult` against time, sampling every `stride` points.
Converts time and energy quantities to the specified units.
Add the new lines to an existing plot and update the legend only if it is the first plot.
If it is not the first plot, add a vertical line to differentiate the segments of the simulation in the plot.
"""
function plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer; first_plot::Bool = false, time_unit = PLOTTING_TIME_UNIT, energy_unit = PLOTTING_ENERGY_UNIT)
    time_range = plotting_time_range(result, stride; time_unit = time_unit)
    if (!first_plot)
        vline!(p, [time_range[1][2]], label = false, color = :black, linestyle = :dot, lw = 2)
    end
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [uconvert(energy_unit, kinetic_energy(result, t)) for (t, time) ∈ time_range],
        label = first_plot ? "Kinetic Energy" : nothing,
        color = 2
    )
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [uconvert(energy_unit, potential_energy(result, t)) for (t, time) ∈ time_range],
        label = first_plot ? "Potential Energy" : nothing,
        color = 1
    )
    plot!(
        p,
        [time for (t, time) ∈ time_range],
        [uconvert(energy_unit, total_energy(result, t)) for (t, time) ∈ time_range],
        label = first_plot ? "Total Energy" : nothing,
        color = 3
    )
end

"""
    plot_rdf(result::MolecularDynamicsResult, σ::Real, start::Integer = 1, stop::Integer = length(result))::Plot

Plot the radial distribution function of a `MolecularDynamicsResult` averaging over the timesteps in `start:stop`
Use `σ` (from Lennard Jones) as a normalization factor for the radius.
"""
function plot_rdf(result::MolecularDynamicsResult, σ::Real, start::Integer = 1, stop::Integer = length(result))
    plot_rdf(result, σ * LENGTH_UNIT, start, stop)
end
"""
    plot_rdf(result::MolecularDynamicsResult, σ::Unitful.Length, start::Integer = 1, stop::Integer = length(result))::Plot

    Plot the radial distribution function of a `MolecularDynamicsResult` averaging over the timesteps in `start:stop`
        Use `σ` (from Lennard Jones) as a normalization factor for the radius.
"""
function plot_rdf(result::MolecularDynamicsResult, σ::Unitful.Length, start::Integer = 1, stop::Integer = length(result))
    N = get_num_bodies(result)
    T = length(result)
    r, g = rdf(result, start, stop)
    @assert length(r) == length(g)
    plot(
        austrip.(r / σ),
        g,
        title = "Radial Distribution Function [n = $(N)] [T = $(T)]",
        xlab = "Distance r/σ",
        ylab = "Radial Distribution g(r)",
        legend = false
    )
end
