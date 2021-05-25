include("../molecular_simulation.jl")
include("../ase_potential_integration.jl")
include("../dftk_integration.jl")
include("../nbs_extensions.jl")

# using Random

include("./nbs_argon.jl")

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

dftk_force_parameters = DFTKForceGenerationParameters(
    box_size,
    ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    [1, 1, 1],
    10u"hartree",
    1e-4
)

# display(analyze_convergence(eq_bodies, force_parameters, [e * u"hartree" for e in (10, 12, 14, 16, 18, 20)]))

# Random.seed!(0)
dftk_forces = generate_forces(eq_bodies, dftk_force_parameters)
println(dftk_forces)

# dftk_force_steps = 100
# dftk_force_stride = dftk_force_steps ÷ 10

# result, bodies = simulate(eq_bodies, dftk_forces, box_size, Δt, dftk_force_steps)
# display(plot_temperature(result, dftk_force_stride))
# display(plot_energy(result, dftk_force_stride))

ase_dftk_force_parameters = ASEForceGenerationParameters(
    box_size,
    ElementCoulomb(:Ar),
    box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    ASEPotential.DFTKCalculatorParameters(
        ecut=10u"hartree",
        kpts=[1, 1, 1],
        scftol=1e-4,
        n_threads=8
    )
)

# Random.seed!(0)
ase_dftk_forces = generate_forces(eq_bodies, ase_dftk_force_parameters)
println(ase_dftk_forces)

;
