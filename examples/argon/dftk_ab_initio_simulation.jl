# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept ab initio MD simulation
# Note that the choice of parameters is for demonstration purposes, and the results are non-physical

using CESMIX
using DFTK
using NBodySimulator
using Unitful
using UnitfulAtomic

N = 8
m = 6.6335209e-26u"kg"
box_size = 1.5u"nm" # this number was chosen arbitrarily
reference_temp = 94.4u"K"
average_v = √(u"k" * reference_temp / m)
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

potential_parameters = LJParameters(
	ϵ = 1.657e-21u"J",
	σ = 0.34u"nm",
	R = 0.765u"nm"
)

initial_bodies = generate_bodies_in_cell_nodes(N, austrip(m), austrip(average_v), austrip(box_size))
eq_parameters = NBSParameters(
	box_size=box_size,
	Δt=Δt,
	steps=20000,
	thermostat=AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
eq_result = @time simulate(initial_bodies, eq_parameters, potential_parameters)

eq_stride = eq_parameters.steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
@time display(plot_rdf(eq_result, potential_parameters.σ, 0.5))

ab_initio_parameters = NBSParameters(
    box_size=box_size,
    Δt=Δt,
    steps=200,
    t₀=eq_parameters.steps * Δt
)
dftk_parameters = DFTKParameters(
    box_size=box_size,
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    Ecut=5u"hartree", # very non-physical but fast for demonstration purposes
    kgrid=[1, 1, 1],
    α=0.7,
    mixing=LdosMixing()
)
ab_initio_result = @time simulate(get_bodies(eq_result), ab_initio_parameters, dftk_parameters)

# Ploting on separate plots because the timespan is so much smaller than in the first phase

ab_initio_stride = ab_initio_parameters.steps ÷ 200

display(plot_temperature(ab_initio_result, ab_initio_stride))
display(plot_energy(ab_initio_result, ab_initio_stride))
@time display(plot_rdf(ab_initio_result, potential_parameters.σ, 1))

write_trajectory(ab_initio_result, box_size, dftk_parameters.psp, dftk_parameters.lattice, "artifacts/argon_ab_initio.traj")

;
