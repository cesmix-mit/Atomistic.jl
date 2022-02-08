var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"This page provides a list of all documented types and functions and in Atomistic.jl.","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"If you are looking for more specifics on the Atomistic API, see the Implementing the Atomistic API page.","category":"page"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [Atomistic]\nOrder   = [:type, :function]","category":"page"},{"location":"api/#Atomistic.MolecularDynamicsResult","page":"API Reference","title":"Atomistic.MolecularDynamicsResult","text":"MolecularDynamicsResult\n\nAbstract type to be extended by all concrete structs representing the result of a molecular dynamics simulation.\n\n\n\n\n\n","category":"type"},{"location":"api/#Atomistic.MolecularDynamicsSimulator","page":"API Reference","title":"Atomistic.MolecularDynamicsSimulator","text":"MolecularDynamicsSimulator\n\nAbstract type to be extended by all concrete molecular dynamics simulators.\n\n\n\n\n\n","category":"type"},{"location":"api/#Atomistic.NBSResult","page":"API Reference","title":"Atomistic.NBSResult","text":"NBSResult <: MolecularDynamicsResult\n\nThe result generated from running an NBSimulator.\n\nField descriptions\n\nresult::SimulationResult the standard simulation result from NBodySimulator\n\n\n\n\n\n","category":"type"},{"location":"api/#Atomistic.NBSimulator","page":"API Reference","title":"Atomistic.NBSimulator","text":"NBSimulator <: MolecularDynamicsSimulator\n\nA wrapper around NBodySimulator to implement the Atomistic API.\n\nField descriptions\n\nΔt::Real the time between timesteps, assumed to be in atomic units\nsteps::Integer the number of timesteps for the simulation\nt₀::Real the starting time of the simulation, assumed to be in atomic units;   defaults to 0\nthermostat::Thermostat the thermostat for the simulation;   many options are defined by NBodySimulator, but a user could also define a custom thermostat;   defaults to the NullThermostat\nsimulator::OrdinaryDiffEqAlgorithm the algorithm to be used for the ODE;   defaults to VelocityVerlet\npotentials::Dict{Symbol,PotentialParameters} dictionary of potentials;   shouldn't be manipulated directly by the user\n\n\n\n\n\n","category":"type"},{"location":"api/#Atomistic.UnimplementedError","page":"API Reference","title":"Atomistic.UnimplementedError","text":"UnimplementedError\n\nException thrown in default implementation of API to indicate that an implementator did not provide an implementation of a particular API function.\n\n\n\n\n\n","category":"type"},{"location":"api/#Atomistic.get_boundary_conditions-Tuple{MolecularDynamicsResult}","page":"API Reference","title":"Atomistic.get_boundary_conditions","text":"get_boundary_conditions(result::MolecularDynamicsResult)::SVector{3, BoundaryCondition}\n\nExtract the boundary conditions of the simulation from the simulation result in an AtomBase-compatible format.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_bounding_box-Tuple{MolecularDynamicsResult}","page":"API Reference","title":"Atomistic.get_bounding_box","text":"get_bounding_box(result::MolecularDynamicsResult)::SVector{3, SVector{3, Unitful.Length}}\n\nExtract the unit-anotated bounding box of the simulation from the simulation result in an AtomBase-compatible format.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_num_bodies-Tuple{MolecularDynamicsResult}","page":"API Reference","title":"Atomistic.get_num_bodies","text":"get_num_bodies(result::MolecularDynamicsResult)::Integer\n\nExtract the number of bodies in the simulation from the simulation result.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_particles-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.get_particles","text":"get_particles(result::MolecularDynamicsResult, t::Integer)::AbstractVector{Atom}\n\nExtract the position of each particle in the system at a particular timestep from the simulation result in an AtomBase-compatible format. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_positions-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.get_positions","text":"get_positions(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Length}}\n\nExtract the unit-anotated position of each particle in the system at a particular timestep from the simulation result in an AtomBase-compatible format. The returned positions should be normalized such that they lie within the bounding box of the system, even if the system is periodic. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_system-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.get_system","text":"get_system(result::MolecularDynamicsResult, t::Integer)::AbstractSystem{3}\n\nExtract the underlying system at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nThe default implementation combines the results of get_particles, get_bounding_box, get_boundary_conditions, and get_time at the timestep to create a FlexibleSystem wrapped in a DyanmicSystem. An implementer of this API could implement a method of this function for their custom result type if it supports a more efficient way to extract the system.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_time-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.get_time","text":"get_time(result::MolecularDynamicsResult, t::Integer)::Unitful.Time\n\nExtract the unit-anotated time of the simulation at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nThe default implementation extracts the time from the result of get_time_range.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_time_range-Tuple{MolecularDynamicsResult}","page":"API Reference","title":"Atomistic.get_time_range","text":"get_time_range(result::MolecularDynamicsResult)::AbstractVector{Unitful.Time}\n\nExtract the unit-anotated time range of the simulation from the simulation result.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.get_velocities-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.get_velocities","text":"get_velocities(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Velocity}}\n\nExtract the unit-anotated velocity of each particle in the system at a particular timestep from the simulation result in an AtomBase-compatible format. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.kinetic_energy-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.kinetic_energy","text":"kinetic_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy\n\nExtract the unit-anotated kinetic energy of the simulation at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.plot_energy!","page":"API Reference","title":"Atomistic.plot_energy!","text":"plot_energy!(p::Plot, result::MolecularDynamicsResult, stride::Integer, first_plot::Bool = false)::Plot\n\nPlot the kinetic, potential, and total energy of a MolecularDynamicsResult against time, sampling every stride points. Add the new lines to an existing plot and update the legend only if it is the first plot. If it is not the first plot, add a vertical line to differentiate the segments of the simulation in the plot.\n\n\n\n\n\n","category":"function"},{"location":"api/#Atomistic.plot_energy-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.plot_energy","text":"plot_energy(result::MolecularDynamicsResult, stride::Integer)::Plot\n\nPlot the kinetic, potential, and total energy of a MolecularDynamicsResult against time, sampling every stride points.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.plot_rdf","page":"API Reference","title":"Atomistic.plot_rdf","text":"plot_rdf(result::MolecularDynamicsResult, σ::Unitful.Length, start::Integer = 1, stop::Integer = length(result))::Plot\n\nPlot the radial distribution function of a `MolecularDynamicsResult` averaging over the timesteps in `start:stop`\n    Use `σ` (from Lennard Jones) as a normalization factor for the radius.\n\n\n\n\n\n","category":"function"},{"location":"api/#Atomistic.plot_rdf-2","page":"API Reference","title":"Atomistic.plot_rdf","text":"plot_rdf(result::MolecularDynamicsResult, σ::Real, start::Integer = 1, stop::Integer = length(result))::Plot\n\nPlot the radial distribution function of a MolecularDynamicsResult averaging over the timesteps in start:stop Use σ (from Lennard Jones) as a normalization factor for the radius.\n\n\n\n\n\n","category":"function"},{"location":"api/#Atomistic.plot_temperature!","page":"API Reference","title":"Atomistic.plot_temperature!","text":"plot_temperature!(p::Plot, result::MolecularDynamicsResult, stride::Integer, first_plot::Bool = false)::Plot\n\nPlot the temperature of a MolecularDynamicsResult against time, sampling every stride points. Adds the new line to an existing plot and update the legend only if it is the first plot. If it is not the first plot, adds a vertical line to differentiate the segments of the simulation in the plot.\n\n\n\n\n\n","category":"function"},{"location":"api/#Atomistic.plot_temperature-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.plot_temperature","text":"plot_temperature(result::MolecularDynamicsResult, stride::Integer)::Plot\n\nPlot the temperature of a MolecularDynamicsResult against time, sampling every stride points.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.potential_energy-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.potential_energy","text":"potential_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy\n\nExtract the unit-anotated potential energy of the simulation at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.rdf-Tuple{MolecularDynamicsResult, Integer, Integer}","page":"API Reference","title":"Atomistic.rdf","text":"rdf(result::MolecularDynamicsResult, start::Integer, stop::Integer)::Tuple{AbstractVector{Unitful.Length},AbstractVector{Real}}\n\nCalculate the radial distribution function from the simulation result.\n\nTo include only a portion of the timesteps for faster computation, increase the start parameter or decrease the stop parameter; by default, the full time range will be iterated over. The result is a named tuple of vectors:     - r: unit-annotated interparticle radial distance bins     - g: distribution value of each bin\n\nThe default implementation provided uses the provided implementations of get_num_bodies, get_bounding_box, get_boundary_conditions, and get_positions. An implementor of the API could use a built in implementation if one that fits this spec is available.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.reference_temperature-Tuple{MolecularDynamicsResult}","page":"API Reference","title":"Atomistic.reference_temperature","text":"reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}\n\nExtract the unit-anotated reference temperature of the simulation from the simulation result. If there is no thermostat with a reference temperature in this simulation, return missing.\n\nAn implementer of this API should implement a method of this function for their custom result type if it supports thermostats. If not implmented, the default implemenation just returns missing.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.simulate-Tuple{AtomsBase.AbstractSystem{3}, MolecularDynamicsSimulator, InteratomicPotentials.ArbitraryPotential}","page":"API Reference","title":"Atomistic.simulate","text":"simulate(system::AbstractSystem{3}, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult\n\nRun a molecular dynamics simulation configured with a particular simulator and potential with any abstract system.\n\nAn implementer of this API should implement a method of this function for their custom simulator type. If the simulator has a fast path for some types of potential, those should be implemented with multiple dispatch.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.temperature-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.temperature","text":"temperature(result::MolecularDynamicsResult, t::Integer)::Unitful.Temperature\n\nExtract the unit-anotated temperature of the simulation at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nAn implementer of this API should implement a method of this function for their custom result type.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.total_energy-Tuple{MolecularDynamicsResult, Integer}","page":"API Reference","title":"Atomistic.total_energy","text":"total_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy\n\nExtract the unit-anotated total energy of the simulation at a particular timestep from the simulation result. Automatically defaults to the end of the simulation when t is not passed.\n\nThe default implementation simply sums the results of kinetic_energy and potential_energy at the timestep. An implementer of this API could implement a method of this function for their custom result type if it supports a more efficient way to calculate this quantity.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.write_ase_trajectory-Tuple{MolecularDynamicsResult, String}","page":"API Reference","title":"Atomistic.write_ase_trajectory","text":"write_ase_trajectory(result::MolecularDynamicsResult, element::DFTK.Element, lattice, filename::String)\n\nWrite the trajectory of a MolecularDynamicsResult to a .traj file.\n\nThe file can be visualized by running ase gui <filename> on the command line.\n\n\n\n\n\n","category":"method"},{"location":"api/#Atomistic.write_nbs_animation-Tuple{NBSResult, String}","page":"API Reference","title":"Atomistic.write_nbs_animation","text":"write_nbs_animation(result::NBSResult, filename::String)\n\nAnimate an NBSResult and store the result in a .gif file.\n\n\n\n\n\n","category":"method"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"EditURL = \"https://github.com/cesmix-mit/Atomistic.jl/blob/main/docs/src/usage.jl\"","category":"page"},{"location":"usage/#Using-Atomistic-Compatible-Packages","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"There are three main steps to using an implementation of the Atomistic API. This simple toy example will walk you through each step using NBodySimulator.jl for the dynamics. Some of the values from this example were taken from this guide.","category":"page"},{"location":"usage/#Step-0:-Load-Dependencies","page":"Using Atomistic-Compatible Packages","title":"Step 0: Load Dependencies","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"using Atomistic\nusing AtomsBase\nusing InteratomicPotentials\nusing NBodySimulator\nusing Unitful\nusing UnitfulAtomic","category":"page"},{"location":"usage/#Step-1A:-Configuring-the-System","page":"Using Atomistic-Compatible Packages","title":"Step 1A: Configuring the System","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"First, we must create an AtomsBase.jl-style system. Here we will create a small Argon cluster with periodic boundary conditions. We could use the AtomsBase constructors directly, but in this case we will use NBodySimulator helper functions to initialize the system.","category":"page"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"N = 864\nelement = :Ar\nbox_size = 3.47786u\"nm\"\nreference_temp = 94.4u\"K\"\n\ninitial_bodies = generate_bodies_in_cell_nodes(N, element, box_size, reference_temp)\ninitial_system = FlexibleSystem(initial_bodies, CubicPeriodicBoundaryConditions(austrip(box_size)))","category":"page"},{"location":"usage/#Step-1B:-Configuring-the-Simulator","page":"Using Atomistic-Compatible Packages","title":"Step 1B: Configuring the Simulator","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"Second, we initialize our simulator. In this case we are using an AndersenThermostat which is provided by NBodySimulator.","category":"page"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"Δt = 1e-2u\"ps\"\nsteps = 2000\nthermostat = AndersenThermostat(austrip(reference_temp), 0.1 / austrip(Δt))\nsimulator = NBSimulator(Δt, steps, thermostat = thermostat)","category":"page"},{"location":"usage/#Step-1C:-Configuring-the-Potential","page":"Using Atomistic-Compatible Packages","title":"Step 1C: Configuring the Potential","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"Lastly, we specify the interatomic potential that we will use for the simulation, Lennard-Jones in this case.","category":"page"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"potential = LennardJones(austrip(1.657e-21u\"J\"), austrip(0.34u\"nm\"), austrip(0.765u\"nm\"), [:Ar])","category":"page"},{"location":"usage/#Step-2:-Running-the-Simulation","page":"Using Atomistic-Compatible Packages","title":"Step 2: Running the Simulation","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"result = simulate(initial_system, simulator, potential)","category":"page"},{"location":"usage/#Step-3:-Analyzing-the-Results","page":"Using Atomistic-Compatible Packages","title":"Step 3: Analyzing the Results","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"We can now analyze the simulation results. The Atomistic API exposes a variety of quantities from each timestep of the simulation (more details here). In this example we will look at the temperature, energy, and radial distribution function.","category":"page"},{"location":"usage/#Temperature","page":"Using Atomistic-Compatible Packages","title":"Temperature","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"plot_temperature(result, simulator.steps ÷ 200)","category":"page"},{"location":"usage/#Energy","page":"Using Atomistic-Compatible Packages","title":"Energy","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"plot_energy(result, simulator.steps ÷ 200)","category":"page"},{"location":"usage/#Radial-Distribution-Function-(RDF)","page":"Using Atomistic-Compatible Packages","title":"Radial Distribution Function (RDF)","text":"","category":"section"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"plot_rdf(result, potential.σ, Int(0.95 * steps))","category":"page"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"","category":"page"},{"location":"usage/","page":"Using Atomistic-Compatible Packages","title":"Using Atomistic-Compatible Packages","text":"This page was generated using Literate.jl.","category":"page"},{"location":"extension/#Implementing-the-Atomistic-API","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"There are two main parts to the Atomistic API. An implementer must provide a concrete simulator type subtyping MolecularDynamicsSimulator and a concrete result type subtyping MolecularDynamicsResult. A summary of the required functionality for each component follows.","category":"page"},{"location":"extension/#Molecular-Dynamics-Simulators","page":"Implementing the Atomistic API","title":"Molecular Dynamics Simulators","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The only function that must be implemented for subtypes of MolecularDynamicsSimulator is a method of simulate.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"simulate(system::AbstractSystem, simulator::MolecularDynamicsSimulator, potential::ArbitraryPotential)::MolecularDynamicsResult","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The first parameter is an instance of an implementation of the AbstractSystem interface defined by AtomsBase.jl. It is important that the implementer rely only on the general API rather than any specific implementation details because the user could provide any system implementation; if there are performance-sensitive decisions based on the system type, this should be achieved through multiple dispatch. The species, position, and velocity of each body as well as the boundary shape and conditions are stored in this data structure. An implementer should take special care to handle units appropriately when converting from this representation to the internal representation of the underlying simulator package.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The second parameter is the custom simulator. This could specify information such as the duration of the simulation, numerical details such as discretization parameters and integration methods, or domain specific options such as thermostats and barostats.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The third parameter is an instance of an implementation of the ArbitraryPotential interface defined by InteratomicPotentials.jl. The main function from this interface that is useful in the dynamics context is InteratomicPotentials.force(s::AbstractSystem, p::ArbitraryPotential)::AbstractVector{StaticVector{3, Real}}. However, it is recommended that implementers consider using InteratomicPotentials.energy_and_force(s::AbstractSystem, p::ArbitraryPotential)::NamedTuple{(:e, :f), Tuple{Real, Vector{SVector{3, Real}}}} and cache the energy results on the MolecularDynamicsResult struct because the additional energy calculation is asymptotically free for most interatomic potential implementations. See the discussion of the potential_energy function below.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The simulate method should return the corresponding implementation of the MolecularDynamicsResult interface as described below.","category":"page"},{"location":"extension/#MolecularDynamicsResult_Specification","page":"Implementing the Atomistic API","title":"Molecular Dynamics Results","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"Fifteen total functions must be implemented for subtypes of MolecularDynamicsResult, four of which have default implementations. The functions allow users to access simulation data from three categories: simulation configuration, time-dependent system data, and simulation analysis.","category":"page"},{"location":"extension/#Simulation-Configuration-Functions","page":"Implementing the Atomistic API","title":"Simulation Configuration Functions","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"There are five functions in the MolecularDynamicsResult API which return time-independent simulation configuration information:","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_time_range(result::MolecularDynamicsResult)::AbstractVector{Unitful.Time}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated vector containing the time value for each step of the simulation which can be iterated for plotting, animation, or other analysis.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_num_bodies(result::MolecularDynamicsResult)::Integer","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the numebr of bodies in the simulation.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_bounding_box(result::MolecularDynamicsResult)::SVector{3, SVector{3, Unitful.Length}}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the AtomsBase.jl-style bounding box of the system used in the simulation.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_boundary_conditions(result::MolecularDynamicsResult)::SVector{3, BoundaryCondition}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the AtomsBase.jl-style boundary conditions of the system used in the simulation.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"reference_temperature(result::MolecularDynamicsResult)::Union{Unitful.Temperature,Missing}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated temperature that describes the temperature maintained by the thermostat used in the simulation. If there is no reference temperature for the particular simulation, the function should return missing, which is the default implementation.","category":"page"},{"location":"extension/#Time-Dependent-System-Data","page":"Implementing the Atomistic API","title":"Time-Dependent System Data","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"There are five functions in the MolecularDynamicsResult API which return system data at a particular timestep, two of which have default implementations.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_time(result::MolecularDynamicsResult, t::Integer)::Unitful.Time","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This convenience function converts a timestep to the unit-annotated time of the timestep using the provided get_time_range implementation. The timestep defaults to the end of the simulation when t is not passed. Most implementers will not need to write a custom implementation.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_positions(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Length}}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the unit-annotated position of each particle in the system at a particular timestep from the simulation result. The returned positions should be normalized such that they lie within the bounding box of the system, even if the system is periodic. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_velocities(result::MolecularDynamicsResult, t::Integer)::AbstractVector{SVector{3, Unitful.Velocity}}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the unit-annotated velocities of each particle in the system at a particular timestep from the simulation result. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_particles(result::MolecularDynamicsResult, t::Integer)::AbstractVector{Atom}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return the particles in the system at a particular timestep from the simulation result. The timestep defaults to the end of the simulation when t is not passed. It is important for any Atom metadata passed in on the input system be preserved when reproducing the particles as some interatomic potential implementations may depend on this data for correctness or performance.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"get_system(result::MolecularDynamicsResult, t::Integer)::AbstractSystem","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return an AtomsBase.jl AbstractSystem which captures the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The default implementation creates a FlexibleSystem by combining the provided implementations of get_particles, get_bounding_box, and get_boundary_conditions.","category":"page"},{"location":"extension/#Simulation-Analysis","page":"Implementing the Atomistic API","title":"Simulation Analysis","text":"","category":"section"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"temperature(result::MolecularDynamicsResult, t::Integer)::Unitful.Temperature","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated temperature for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"kinetic_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated kinetic energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"potential_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated potential energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when t is not passed. If potential energy values were not cached at simulation time, the relevant function from the InteratomicPotentials.jl interface is InteratomicPotentials.potential_energy(a::AbstractSystem, p::ArbitraryPotential)::Real.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"total_energy(result::MolecularDynamicsResult, t::Integer)::Unitful.Energy","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should return a unit-anotated total energy for the system at a particular timestep in the simulation. The timestep defaults to the end of the simulation when t is not passed.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The default implementation simply sums the kinetic and potential energy functions, but an implementor of the API might provide a custom implementation if there is a more direct means of calculation provided by the underlying simulator.","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"rdf(result::MolecularDynamicsResult, start::Integer, stop::Integer)::Tuple{AbstractVector{Real},AbstractVector{Real}}","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"This function should calculate the radial distribution function of the system averaged across a range of timesteps. It should return a named tuple of vectors:","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"r: unit-annotated interparticle radial distance bins\ng: distribution value of each bin","category":"page"},{"location":"extension/","page":"Implementing the Atomistic API","title":"Implementing the Atomistic API","text":"The default implementation provided uses the provided implementations of get_num_bodies, get_bounding_box, get_boundary_conditions, and get_positions. An implementor of the API could use a built in implementation if one that fits this spec is available.","category":"page"},{"location":"#[WIP]-Atomistic.jl","page":"Home","title":"[WIP] Atomistic.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Package that provides an integrated workflow for molecular dynamics (MD) simulations. Defines an API for MD simulations that is compatible with the interatomic potential interface defined by InteratomicPotentials.jl and the atomic configuration interface defined by AtomsBase.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Developed as part of the CESMIX Julia package suite. See also InteratomicPotentials.jl, PotentialLearning.jl, and PotentialUQ.jl.","category":"page"},{"location":"#Conventions","page":"Home","title":"Conventions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The unit convention throughout the package and other packages in the CESMIX Julia Suite is to assume all unspecified units to be atomic units as defined in the UnitfulAtomic.jl package. All exposed interfaces should allow for numeric or unitful input. For clarity's sake, it is strongly recommended that user code utilize Unitful wherever possible. Internally, Atomistic will automatically convert these quantities to be compatible without a significant performance penalty.","category":"page"},{"location":"#Next-Steps","page":"Home","title":"Next Steps","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you want to integrate an existing MD code with the Atomistic API, see Implementing the Atomistic API. If you want to use a code that is already integrated with Atomistic to run MD simulations, see Using Atomistic-Compatible Packages. If you want to see the full API reference (for Atomistic API and other types and functions exported by the package), see API Reference.","category":"page"}]
}
