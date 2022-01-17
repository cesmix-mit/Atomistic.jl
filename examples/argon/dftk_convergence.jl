# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

using Atomistic
using AtomsBase
using DFTK
using NBodySimulator
using Unitful
using UnitfulAtomic

setup_threading(n_blas = 4)

N = 8
element = :Ar
box_size = 1.5u"nm" # this number was chosen arbitrarily
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_bodies = generate_bodies_in_cell_nodes(N, element, box_size, reference_temp)
initial_system = FlexibleSystem(initial_bodies, CubicPeriodicBoundaryConditions(austrip(box_size)))

eq_steps = 20000
eq_thermostat = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
eq_simulator = NBSimulator(Δt, eq_steps, thermostat = eq_thermostat)
potential = LennardJonesParameters(1.657e-21u"J", 0.34u"nm", 0.765u"nm")

eq_result = @time simulate(initial_system, eq_simulator, potential)

display(@time plot_rdf(eq_result, potential.σ, 0.5))

dftk_potential = DFTKPotential(
    psp = ElementPsp(:Ar, psp = load_psp(list_psp(:Ar, functional = "lda")[1].identifier)),
    lattice = box_size * [[1.0 0 0]; [0 1.0 0]; [0 0 1.0]],
    Ecut = 5u"hartree",
    kgrid = [1, 1, 1],
    damping = 0.7,
    mixing = LdosMixing()
)

display(@time analyze_convergence(get_system(eq_result), dftk_potential, (5:2.5:25)u"hartree"))

;
