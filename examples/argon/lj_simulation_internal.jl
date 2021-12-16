# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf

using Atomistic
using AtomsBase
using NBodySimulator
using Plots
using Unitful
using UnitfulAtomic

N = 864
element = :Ar
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_bodies = generate_bodies_in_cell_nodes(N, element, box_size, reference_temp)
initial_system = DynamicSystem(initial_bodies, box_size)
eq_simulator = NBSimulator(
    Δt = Δt,
    steps = 2000,
    t₀ = 0.0u"s",
    thermostat = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
potential = LJPotential(
    ϵ = 1.657e-21u"J",
    σ = 0.34u"nm",
    R = 0.765u"nm"
)
eq_result = @time simulate(initial_system, eq_simulator, potential)

temp = @time plot_temperature(eq_result, 10)
energy = @time plot_energy(eq_result, 10)

prod_simulator = NBSimulator(
    Δt = Δt,
    steps = 5000,
    t₀ = eq_simulator.steps * Δt
)
prod_result = @time simulate(get_system(eq_result), prod_simulator, potential)

display(@time plot_temperature!(temp, prod_result, 10))
display(@time plot_energy!(energy, prod_result, 10))

rdf = @time plot_rdf(prod_result, potential.σ, 0.05)
display(rdf)
savefig(rdf, "artifacts/argon_lj_rdf.svg")

;
