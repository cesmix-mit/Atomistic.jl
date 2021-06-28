# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

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
	thermostat=argon_equilibration_thermostat(reference_temp, thermostat_prob)
)
eq_result, eq_bodies = simulate(initial_bodies, eq_parameters)

display(plot_rdf(eq_result, sample_fraction=2))

dftk_parameters = DFTKParameters(
    box_size=box_size,
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    Ecut=10u"hartree",
    kgrid=[1, 1, 1],
    α=0.7,
    mixing=LdosMixing()
)

display(analyze_convergence(eq_bodies, dftk_parameters, [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25]u"hartree"))

;
