m = 6.6335209e-26u"kg"
σ = 0.34u"nm"
cutoff = 0.765u"nm"
ϵ = 1.657e-21u"J"

argon_initial_bodies(N::Integer, box_size::Quantity, reference_temp::Quantity) = generate_bodies_in_cell_nodes(N, austrip(m), austrip(√(u"k" * reference_temp / m)), austrip(box_size))
argon_lennard_jones() = Dict(:lennard_jones => LennardJonesParameters(austrip(ϵ), austrip(σ), austrip(cutoff)))
argon_equilibration_thermostat(reference_temp::Quantity, thermostat_prob::AbstractFloat, Δt::Quantity) = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
