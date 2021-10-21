# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

using Atomistic
using InteratomicPotentials
using NBodySimulator
using Plots
using Unitful
using UnitfulAtomic

N = 864
m = 6.6335209e-26u"kg"
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"
average_v = √(u"k" * reference_temp / m)
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

potential = LennardJones(
	austrip(1.657e-21u"J"),
	austrip(0.34u"nm")
)

initial_bodies = MassBodies(N, m, average_v, box_size)
eq_simulator = NBSimulator(
	Δt=Δt,
	steps=200,
	t₀=0.0u"s",
	thermostat=AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
eq_result = @time simulate(initial_bodies, eq_simulator, potential)

eq_stride = 10

temp = @time plot_temperature(eq_result, eq_stride)
energy = @time plot_energy(eq_result, eq_stride)

prod_simulator = NBSimulator(
	Δt=Δt,
	steps=500,
	t₀=eq_simulator.steps * Δt
)
prod_result = @time simulate(get_bodies(eq_result), prod_simulator, potential)

prod_stride = 10

@time plot_temperature!(temp, prod_result, prod_stride)
display(temp)
@time plot_energy!(energy, prod_result, prod_stride)
display(energy)

# rdf = @time plot_rdf(prod_result, potential.σ, 0.05)
# display(rdf)
# savefig(rdf, "artifacts/argon_lj_rdf.svg")

;
