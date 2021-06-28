# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

include("../../src/molecular_simulation.jl")

include("./nbs_argon.jl")

N = 864
box_size = 3.47786u"nm"

reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily

eq_steps = 2000
Δt = 1e-2u"ps"

eq_result, eq_bodies = simulate_lennard_jones_argon_equilibration(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

eq_stride = 10

temp = plot_temperature(eq_result, eq_stride)
energy = plot_energy(eq_result, eq_stride)

prod_steps = 5000

prod_result, prod_bodies = simulate_lennard_jones_argon_production(eq_bodies, box_size, Δt, prod_steps, eq_steps * Δt)

prod_stride = 10

display(plot_temperature!(temp, prod_result, prod_stride))
display(plot_energy!(energy, prod_result, prod_stride))

rdf = plot_rdf(prod_result)
savefig(rdf, "artifacts/argon_lj_rdf.svg")
display(rdf)

;
