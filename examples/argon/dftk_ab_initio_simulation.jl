# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept ab initio MD simulation
# Note that the choice of parameters is for demonstration purposes and the results are non-physical

include("../../src/molecular_simulation.jl")

include("./nbs_argon.jl")

N = 8
box_size = 4σ # arbitrarly choosing 4σ
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_bodies = argon_initial_bodies(N, box_size, reference_temp)
eq_parameters = NBSParameters(
	potentials=argon_lennard_jones(),
	box_size=box_size,
	Δt=Δt,
	steps=20000,
	thermostat=argon_equilibration_thermostat(reference_temp, thermostat_prob, Δt)
)
eq_result, eq_bodies = simulate(initial_bodies, eq_parameters)

eq_stride = eq_parameters.steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
display(plot_rdf(eq_result, sample_fraction=2))

nbs_parameters = NBSParameters(
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
ab_initio_result, ab_initio_bodies = simulate(eq_bodies, nbs_parameters, dftk_parameters)

# Ploting on separate plots because the timespan is so much smaller than in the first phase

ab_initio_stride = 1

display(plot_temperature(ab_initio_result, ab_initio_stride))
display(plot_energy(ab_initio_result, ab_initio_stride))
display(plot_rdf(ab_initio_result, σ=σ, sample_fraction=1))

write_trajectory(ab_initio_result, box_size, dftk_parameters.psp, dftk_parameters.lattice, "artifacts/argon_ab_initio.traj")

;
