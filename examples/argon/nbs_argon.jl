function argon_simulate_equilibration(N::Integer, box_size::Quantity, Δt::Quantity, steps::Integer, reference_temp::Quantity, thermostat_prob::AbstractFloat)
    m = auconvert(6.6335209e-26u"kg")
    mean_v = auconvert(√(u"k" * reference_temp / m))

    initial_bodies = generate_bodies_in_cell_nodes(N, ustrip(m), ustrip(mean_v), ustrip(box_size))
	thermostat = AndersenThermostat(ustrip(reference_temp), thermostat_prob / ustrip(Δt))

	return simulate(initial_bodies, argon_lennard_jones(), box_size, Δt, steps, 0.0u"ps", thermostat)
end

function argon_simulate_production(bodies::Vector{MassBody}, box_size::Quantity, Δt::Quantity, steps::Integer, t_0::Quantity)
	return simulate(bodies, argon_lennard_jones(), box_size, Δt, steps, t_0)
end

function argon_lennard_jones()
	σ = auconvert(0.34u"nm")
	cutoff = auconvert(0.765u"nm")
	ϵ = auconvert(1.657e-21u"J")
	return Dict(:lennard_jones => LennardJonesParameters(ustrip(ϵ), ustrip(σ), ustrip(cutoff)))
end