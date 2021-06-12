# Script to demonstrate convergence analysis of ecut values for DFTK.jl

include("../../src/molecular_simulation.jl")
include("../../src/dftk_integration.jl")
include("../../src/nbs_extensions.jl")

include("./nbs_argon.jl")

N = 5
box_size = auconvert(3.47786u"nm")

reference_temp = auconvert(94.4u"K")
thermostat_prob = 0.1

eq_steps = 100000
Δt = auconvert(1e-2u"ps")

eq_result, eq_bodies = argon_simulate_equilibration(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

display(plot_rdf(eq_result))

dftk_force_parameters = DFTKForceGenerationParameters(
    box_size,
    ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    [1, 1, 1],
    10u"hartree",
    1e-4
)

display(analyze_convergence(eq_bodies, dftk_force_parameters, [e * u"hartree" for e in (10, 12, 14, 16, 18, 20)]))

;
