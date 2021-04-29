include("./nbs_argon.jl")
include("./dftk_argon.jl")

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

# display(analyze_convergence(eq_bodies, box_size))

forces = calculate_dftk_forces(eq_bodies, box_size)

dftk_force_steps = 100
dftk_force_stride = dftk_force_steps ÷ 10

result, bodies = dftk_force_simulate(eq_bodies, forces, box_size, Δt, dftk_force_steps)
display(plot_temperature(result, dftk_force_stride))
display(plot_energy(result, dftk_force_stride))

;
