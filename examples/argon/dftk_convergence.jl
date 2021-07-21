# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

using CESMIX
using DFTK
using NBodySimulator
using Unitful
using UnitfulAtomic

N = 8
m = 6.6335209e-26u"kg"
box_size = 1.5u"nm" # this number was chosen arbitrarily
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

potential_parameters = LJParameters(
	ϵ = 1.657e-21u"J",
	σ = 0.34u"nm",
	R = 0.765u"nm"
)

initial_bodies = generate_bodies_in_cell_nodes(N, austrip(m), austrip(√(u"k" * reference_temp / m)), austrip(box_size))
eq_parameters = NBSParameters(
	box_size=box_size,
	Δt=Δt,
	steps=20000,
	thermostat=AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
)
eq_result, eq_bodies = simulate(initial_bodies, eq_parameters, potential_parameters)

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
