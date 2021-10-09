# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

using Atomistic
using DFTK
using NBodySimulator
using Unitful
using UnitfulAtomic

# setup_threading(n_blas=4)

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

@time display(plot_rdf(eq_result, potential.σ, 0.5))

dftk_potential = DFTKPotential(
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    Ecut=10u"hartree",
    kgrid=[1, 1, 1],
    damping=0.7,
    mixing=LdosMixing()
)

@time display(analyze_convergence(get_bodies(eq_result), dftk_potential, (5:2.5:25)u"hartree"))

;
