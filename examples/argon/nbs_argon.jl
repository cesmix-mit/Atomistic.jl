m = 6.6335209e-26u"kg"
σ = 0.34u"nm"
cutoff = 0.765u"nm"
ϵ = 1.657e-21u"J"

function simulate_lennard_jones_argon_equilibration(N::Integer, box_size::Quantity, Δt::Quantity, steps::Integer, reference_temp::Quantity, thermostat_prob::AbstractFloat)
    mean_v = √(u"k" * reference_temp / m)

    initial_bodies = generate_bodies_in_cell_nodes(N, austrip(m), austrip(mean_v), austrip(box_size))
	thermostat = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))

	return simulate(initial_bodies, lennard_jones_argon(), box_size, Δt, steps, 0.0u"s", thermostat)
end

function simulate_lennard_jones_argon_production(bodies::Vector{MassBody}, box_size::Quantity, Δt::Quantity, steps::Integer, t_0::Quantity)
	return simulate(bodies, lennard_jones_argon(), box_size, Δt, steps, t_0)
end

lennard_jones_argon() = Dict(:lennard_jones => LennardJonesParameters(austrip(ϵ), austrip(σ), austrip(cutoff)))

;