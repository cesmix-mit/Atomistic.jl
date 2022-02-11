# Script to demonstrate convergence analysis of Ecut values for DFTK.jl

using Atomistic
using DFTK
using InteratomicPotentials
using NBodySimulator
using Plots
using Unitful
using UnitfulAtomic

setup_threading(n_blas = 4)

N = 8
element = :Ar
box_size = 1.5u"nm" # this number was chosen arbitrarily
reference_temp = 94.4u"K"
thermostat_prob = 0.1 # this number was chosen arbitrarily
Δt = 1e-2u"ps"

initial_system = generate_atoms_in_cubic_cell(N, element, box_size, reference_temp)
pspkey = list_psp(:Ar, functional = "lda")[1].identifier
for atom ∈ initial_system
    atom.data[:pseudopotential] = pspkey
end

eq_steps = 20000
eq_thermostat = AndersenThermostat(austrip(reference_temp), thermostat_prob / austrip(Δt))
eq_simulator = NBSimulator(Δt, eq_steps, thermostat = eq_thermostat)
potential = LennardJonesParameters(1.657e-21u"J", 0.34u"nm", 0.765u"nm")

eq_result = @time simulate(initial_system, eq_simulator, potential)

display(@time plot_rdf(eq_result, potential.σ))

system = get_system(eq_result)
cutoffs = (5:2.5:25)u"hartree"
energies = @time map(cutoffs) do cutoff
    @info "Ecut: $(cutoff)"
    dftk_potential = DFTKPotential(cutoff, [1, 1, 1]; damping = 0.7)
    InteratomicPotentials.potential_energy(system, dftk_potential)
end

display(plot(
    title = "DFTK Convergence Analysis",
    xlab = "Ecut",
    ylab = "Total Energy",
    legend = false,
    cutoffs,
    energies * u"hartree"
))

;
