# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

include("../../src/molecular_simulation.jl")

include("./nbs_argon.jl")

N = 8
box_size = 4σ # arbitrarly choosing 4σ

reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily

eq_steps = 20000
Δt = 1e-2u"ps"

eq_result, eq_bodies = simulate_lennard_jones_argon_equilibration(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

display(plot_rdf(eq_result, sample_fraction=2))

dftk_parameters = DFTKForceGenerationParameters(
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
