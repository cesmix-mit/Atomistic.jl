# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

using CESMIX
using NBodySimulator
using Plots
using Unitful
using UnitfulAtomic

N = 864
m = 6.6335209e-26u"kg"
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

potential_parameters = LJParameters(
	ϵ = 1.657e-21u"J",
	σ = 0.34u"nm",
	R = 0.765u"nm"
)

initial_bodies = generate_bodies_in_cell_nodes(N, austrip(m), austrip(√(u"k" * reference_temp / m)), austrip(box_size))
eq_parameters = NBSParameters(
	box_size=box_size,
	Δt=Δt,
	steps=2000,
	thermostat=AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
eq_result, eq_bodies = simulate(initial_bodies, eq_parameters, potential_parameters)

eq_stride = 10

temp = plot_temperature(eq_result, eq_stride)
energy = plot_energy(eq_result, eq_stride)

prod_parameters = NBSParameters(
	box_size=box_size,
	Δt=Δt,
	steps=5000,
	t₀=eq_parameters.steps * Δt
)
prod_result, prod_bodies = simulate(eq_bodies, prod_parameters, potential_parameters)

prod_stride = 10

display(plot_temperature!(temp, prod_result, prod_stride))
display(plot_energy!(energy, prod_result, prod_stride))

rdf = plot_rdf(prod_result)
display(rdf)
savefig(rdf, "artifacts/argon_lj_rdf.svg")

;
