# Script to compare the results of interfacing with DFTK.jl directly and through the ASEPotential.jl --> ase.py --> asedftk.py --> DFTK.jl pipeline

using ASEPotential
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

lattice = box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]]
kpts = [1, 1, 1]
ecut = 10u"hartree"
tol=1e-6

dftk_parameters = DFTKParameters(
    box_size=box_size,
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=lattice,
    Ecut=ecut,
    kgrid=kpts,
    tol=tol,
    # α=0.7,
    # mixing=LdosMixing()
)

dftk_forces = generate_forces(eq_bodies, dftk_parameters)
println(dftk_forces)

ase_dftk_parameters = ASEPotentialParameters(
    box_size,
    ElementCoulomb(:Ar),
    lattice,
    ASEPotential.DFTKCalculatorParameters(
        ecut=ecut,
        kpts=kpts,
        scftol=tol,
        n_threads=8
    )
)

ase_dftk_forces = generate_forces(eq_bodies, ase_dftk_parameters)
println(ase_dftk_forces)

# Comparing these results isn't really meaningful because the random number generators aren't seeded in the same way

;
