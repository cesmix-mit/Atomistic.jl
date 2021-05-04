include("./nbs_argon.jl")
include("../dftk_integration.jl")

# N = 864
N = 5
box_size = auconvert(3.47786u"nm")

reference_temp = auconvert(94.4u"K")
thermostat_prob = 0.1

# eq_steps = 2000
eq_steps = 100000
Δt = auconvert(1e-2u"ps")

eq_stride = eq_steps ÷ 200

eq_result, eq_bodies = equilibrate(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
display(plot_rdf(eq_result))

force_parameters = DFTKForceGenerationParameters(
    box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    [1, 1, 1],
    10u"hartree",
    1e-4
)

# display(analyze_convergence(eq_bodies, box_size, force_parameters, [e * u"hartree" for e in (10, 12, 14, 16, 18, 20)]))

forces = generate_forces(eq_bodies, box_size, force_parameters)

dftk_force_steps = 100
dftk_force_stride = dftk_force_steps ÷ 10

result, bodies = simulate(eq_bodies, forces, box_size, Δt, dftk_force_steps)
display(plot_temperature(result, dftk_force_stride))
display(plot_energy(result, dftk_force_stride))

;
