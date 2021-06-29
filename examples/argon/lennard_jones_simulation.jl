# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

include("./nbs_argon.jl")

N = 864
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_bodies = argon_initial_bodies(N, box_size, reference_temp)
eq_parameters = NBSParameters(
	potentials=argon_lennard_jones(),
	box_size=box_size,
	Δt=Δt,
	steps=2000,
	thermostat=argon_equilibration_thermostat(reference_temp, thermostat_prob, Δt)
)
eq_result, eq_bodies = simulate(initial_bodies, eq_parameters)

eq_stride = 10

temp = plot_temperature(eq_result, eq_stride)
energy = plot_energy(eq_result, eq_stride)

prod_parameters = NBSParameters(
	potentials=argon_lennard_jones(),
	box_size=box_size,
	Δt=Δt,
	steps=5000,
	t₀=eq_parameters.steps * Δt
)
prod_result, prod_bodies = simulate(eq_bodies, prod_parameters)

prod_stride = 10

display(plot_temperature!(temp, prod_result, prod_stride))
display(plot_energy!(energy, prod_result, prod_stride))

rdf = plot_rdf(prod_result)
display(rdf)
savefig(rdf, "artifacts/argon_lj_rdf.svg")

;
