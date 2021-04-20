include("./nbs_argon.jl")
# include("./dftk_argon.jl")

# N = 864
N = 20
box_size = auconvert(3.47786u"nm")

reference_temp = auconvert(94.4u"K")
thermostat_prob = 0.1

steps = 2000
# steps = 100000
Δt = auconvert(1e-2u"ps")

stride = steps ÷ 200

result, bodies = equilibrate(N, box_size, Δt, steps, reference_temp, thermostat_prob)

display(plot_temperature(result, stride))
display(plot_energy(result, stride))
display(plot_rdf(result))

# calculate_forces(bodies, box_size)