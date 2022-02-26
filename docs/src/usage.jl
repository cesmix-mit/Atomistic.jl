# # Using Atomistic-Compatible Packages

# There are three main steps to using an implementation of the Atomistic API.
# This simple toy example will walk you through each step using [NBodySimulator.jl](https://github.com/SciML/NBodySimulator.jl) for the dynamics.
# Some of the values from this example were taken from [this guide](https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf).

# ## Step 0: Load Dependencies

using Atomistic
using InteratomicPotentials
using Molly

# ## Step 1A: Configuring the System
# First, we must create an [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl)-style system.
# Here we will create a small Argon cluster with periodic boundary conditions.
# We are using a helper function to initialize the system automatically with a particular temperature.

N = 864
element = :Ar
box_size = 3.47786u"nm"
reference_temp = 94.4u"K"

initial_system = generate_atoms_in_cubic_cell(N, element, box_size, reference_temp)

# ## Step 1B: Configuring the Simulator
# Second, we initialize our simulator. In this case we are using an AndersenThermostat which is provided by NBodySimulator.

Δt = 1e-2u"ps"
steps = 2000
thermostat = Molly.AndersenThermostat(reference_temp, Δt * 10)
simulator = MollySimulator(Δt, steps, coupling = thermostat)

# ## Step 1C: Configuring the Potential
# Lastly, we specify the interatomic potential that we will use for the simulation, Lennard-Jones in this case.

potential = InteratomicPotentials.LennardJones(austrip(1.657e-21u"J"), austrip(0.34u"nm"), austrip(0.765u"nm"), [:Ar])

# ## Step 2: Running the Simulation

result = simulate(initial_system, simulator, potential)

# ## Step 3: Analyzing the Results
# We can now analyze the simulation results.
# The Atomistic API exposes a variety of quantities from each timestep of the simulation (more details [here](@ref MolecularDynamicsResult_Specification)).
# In this example we will look at the temperature, energy, and radial distribution function.

# ### Temperature
plot_temperature(result, simulator.steps ÷ 200)

# ### Energy
plot_energy(result, simulator.steps ÷ 200)

# ### Radial Distribution Function (RDF)
plot_rdf(result, potential.σ, Int(0.95 * steps))
