# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept ab initio MD simulation
# Note that the choice of parameters is for demonstration purposes, and the results are non-physical

using Atomistic
using DFTK
using NBodySimulator
using Unitful
using UnitfulAtomic

setup_threading(n_blas=4)

N = 8
m = 6.6335209e-26u"kg"
box_size = 1.5u"nm" # this number was chosen arbitrarily
reference_temp = 94.4u"K"
average_v = √(u"k" * reference_temp / m)
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

potential = LJParameters(
	ϵ=1.657e-21u"J",
	σ=0.34u"nm",
	R=0.765u"nm"
)

initial_bodies = MassBodies(N, m, average_v, box_size)
eq_simulator = NBSimulator(
	Δt=Δt,
	steps=20000,
	thermostat=AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
eq_result = @time simulate(initial_bodies, eq_simulator, potential)

eq_stride = eq_simulator.steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
@time display(plot_rdf(eq_result, potential.σ, 0.5))

ab_initio_simulator = NBSimulator(
    Δt=Δt,
    steps=200,
    t₀=eq_simulator.steps * Δt
)
dftk_potential = DFTKPotential(
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    Ecut=5u"hartree", # very non-physical but fast for demonstration purposes
    kgrid=[1, 1, 1],
    damping=0.7,
    mixing=LdosMixing()
)
ab_initio_result = @time simulate(get_bodies(eq_result), ab_initio_simulator, dftk_potential)

# Ploting on separate plots because the timespan is so much smaller than in the first phase

ab_initio_stride = ab_initio_simulator.steps ÷ 200

display(plot_temperature(ab_initio_result, ab_initio_stride))
display(plot_energy(ab_initio_result, ab_initio_stride))
@time display(plot_rdf(ab_initio_result, potential.σ, 1))

write_nbs_animation(ab_initio_result, "artifacts/argon_ab_initio.gif")
write_ase_trajectory(ab_initio_result, dftk_potential.psp, dftk_potential.lattice, "artifacts/argon_ab_initio.traj")

;
