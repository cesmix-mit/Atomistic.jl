# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept
# Note that the choice of parameters is for demonstration purposes and the results are non-physical

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

eq_stride = eq_steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
display(plot_rdf(eq_result))

# This scipt only runs DFTK once as a proof of concept -- note that this is not sensible for an ab initio simulation

dftk_force_parameters = DFTKForceGenerationParameters(
    box_size,
    ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    [1, 1, 1],
    10u"hartree",
    1e-4
)

dftk_forces = generate_forces(eq_bodies, dftk_force_parameters)

dftk_force_steps = 100

result, bodies = simulate(eq_bodies, dftk_forces, box_size, Δt, dftk_force_steps)

dftk_force_stride = dftk_force_steps ÷ 10

# Ploting on separate plots because the timespan is so much smaller than in the first phase

display(plot_temperature(result, dftk_force_stride))
display(plot_energy(result, dftk_force_stride))

;
