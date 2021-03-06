# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

using Atomistic
using InteratomicPotentials
using Molly
using Plots

N = 864
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"
coupling_factor = 10 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_system = generate_atoms_in_cubic_cell(N, :Ar, box_size, reference_temp)

eq_steps = 2000
eq_thermostat = Molly.AndersenThermostat(reference_temp, Δt * coupling_factor)
eq_simulator = MollySimulator(Δt, eq_steps, coupling=eq_thermostat)
potential = InteratomicPotentials.LennardJones(1.657e-21u"J", 0.34u"nm", 0.765u"nm", [:Ar])

eq_result = @time simulate(initial_system, eq_simulator, potential)

temp = plot_temperature(eq_result, 10)
energy = plot_energy(eq_result, 10)

prod_steps = 5000
prod_simulator = MollySimulator(Δt, prod_steps, t₀=get_time(eq_result))

prod_result = @time simulate(get_system(eq_result), prod_simulator, potential)

display(plot_temperature!(temp, prod_result, 10))
display(plot_energy!(energy, prod_result, 10))

rdf = plot_rdf(prod_result, potential.σ, Int(0.95 * prod_steps))
display(rdf)
savefig(rdf, "artifacts/argon_lj_molly_ip_rdf.svg")

;
