# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept ab initio MD simulation
# Note that the choice of parameters is for demonstration purposes and the results are non-physical

include("../../src/molecular_simulation.jl")

include("./nbs_argon.jl")

N = 8
box_size = 4σ # arbitrarly choosing 4σ

reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily

eq_steps = 20000
Δt = 1e-2u"ps"

eq_result, eq_bodies = simulate_lennard_jones_argon_equilibration(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

eq_stride = eq_steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
display(plot_rdf(eq_result, sample_fraction=2))

ab_initio_parameters = AbInitioPotentialParameters(
    forceGenerationParameters=DFTKForceGenerationParameters(
        box_size=box_size,
        psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
        lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
        Ecut=5u"hartree", # very non-physical but fast for demonstration purposes
        kgrid=[1, 1, 1],
        α=0.7,
        mixing=LdosMixing()
    )
)

ab_initio_steps = 200

ab_initio_result, ab_initio_bodies = simulate(eq_bodies, ab_initio_parameters, box_size, Δt, ab_initio_steps)

ab_initio_stride = 1

# Ploting on separate plots because the timespan is so much smaller than in the first phase

display(plot_temperature(ab_initio_result, ab_initio_stride))
display(plot_energy(ab_initio_result, ab_initio_stride))
display(plot_rdf(ab_initio_result, σ=σ, sample_fraction=1))

write_trajectory(
    ab_initio_result,
    box_size,
    ab_initio_parameters.forceGenerationParameters.psp,
    ab_initio_parameters.forceGenerationParameters.lattice,
    "artifacts/argon_ab_initio.traj"
)

;