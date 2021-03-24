### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ b148ddb6-817b-11eb-2878-c1b583bb35f2
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			"LinearAlgebra"
			"Molly"
			"Parameters"
			"Plots"
	])
	using LinearAlgebra
	using Molly
	using Plots
end

# ╔═╡ 3328ea8e-8cbb-11eb-2c3f-b143741525be
# source: https://github.com/fonsp/Pluto.jl/issues/115
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end;

# ╔═╡ 4ee38c78-8b4b-11eb-3284-45e3dd145d1f
begin
	M = ingredients("../src/molly_extensions.jl");
	import .M:
		FixedVelocityThermostat,
		KineticEnergyLogger,
		PotentialEnergyLogger,
		VelocityLogger
end

# ╔═╡ d8daea28-817a-11eb-26e1-a909bc31bf83
md"""
# Example MD Simulation of Argon Gas
"""

# ╔═╡ c7559cde-817b-11eb-1d86-b5df501583b5
md"""
This is a basic simulation of a cluster of Argon molecules in a cubical box. The purpose of this notebook is to experiment with Molly.jl with the built-in Lennard-Jones potential and demonstrate how to use the package to reproduce a basic known experiment.
"""

# ╔═╡ baaabfde-817c-11eb-2c74-bbcff5f1d162
md"""
## Parameters
"""

# ╔═╡ cd5a9930-817c-11eb-0e02-cd5a7462b687
md"""
The parameters for this example are from [this guide](https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf).
"""

# ╔═╡ d710689c-817c-11eb-0727-3b104fbbb669
begin
	n_atoms = 864
	mass = 6.6335209e-26 / 1.66054e-27 # amu
	
	box_size = 3.47786 # nm
	
	σ = 0.34 # nm
	cutoff = 0.765 # nm
	ϵ = 1.657e-21 / 1.602176634e-19 # eV
	
	reference_temp = 94.4 # K
	
	eq_steps = 2000
	prod_steps = 5000
	timestep = 1e-2 # ps
	stride = 10
end;

# ╔═╡ 04a5c10c-8cca-11eb-1386-b128e8ce0ed8
md"""
## Equilibration Stage
"""

# ╔═╡ b4379f32-817e-11eb-2689-9d11da52ae75
md"""
### Initialization
"""

# ╔═╡ cfc54920-817e-11eb-2758-4748397ddf98
md"""
The referenced guide uses a GROMACS coordinate file for the initialization, but we are just using a random distribution over the volume of the box.
"""

# ╔═╡ 920c3baa-817e-11eb-2037-410965f2894d
begin
	simulator = VelocityVerlet()
	atoms = [Atom(mass=mass, σ=σ, ϵ=ϵ) for i ∈ 1:n_atoms]
	general_inters = (LennardJones(ShiftedPotentialCutoff(cutoff), true),)
	# general_inters = (LennardJones(),)
	initial_coords = [box_size .* rand(SVector{3}) for i ∈ 1:n_atoms]
	initial_velocities = [velocity(mass, reference_temp) for i ∈ 1:n_atoms]
	neighbour_finder = DistanceNeighbourFinder(trues(n_atoms, n_atoms), 1, cutoff)
	eq_thermostat = FixedVelocityThermostat()
	# eq_thermostat = AndersenThermostat(1.0)
	eq_loggers = Dict(
		"temp" => TemperatureLogger(stride),
		"coords" => CoordinateLogger(stride),
		"energy" => EnergyLogger(stride),
		"kinetic" => KineticEnergyLogger(stride),
		"potential" => PotentialEnergyLogger(stride),
		"velocities" => VelocityLogger(stride)
	)
	eq_simulation = Simulation(
		simulator=simulator,
		atoms=atoms,
		general_inters=general_inters,
		coords=deepcopy(initial_coords),
		velocities=deepcopy(initial_velocities),
		temperature=reference_temp,
		box_size=box_size,
		neighbour_finder=neighbour_finder,
		thermostat=eq_thermostat,
		loggers=eq_loggers,
		timestep=timestep,
		n_steps=eq_steps
	);
end;

# ╔═╡ 1c301082-8180-11eb-158f-7538d52db5c3
md"""
### Equilibration Simulation
"""

# ╔═╡ d2fa5508-817f-11eb-1f96-19b27221e7be
begin
	simulate!(eq_simulation, parallel=true)
	eq_time_range = collect(1:eq_steps/stride) .* timestep
	eq_temperatures = eq_simulation.loggers["temp"].temperatures
	eq_energies = eq_simulation.loggers["energy"].energies
	eq_potential_energies = eq_simulation.loggers["potential"].energies
	eq_kinetic_energies = eq_simulation.loggers["kinetic"].energies
	eq_coords = eq_simulation.loggers["coords"].coords[end]
	eq_velocities = eq_simulation.loggers["velocities"].velocities[end]
end;

# ╔═╡ 535b18c4-8cc3-11eb-3c09-c5118d0a2d2f
norm.(initial_velocities)

# ╔═╡ 5a247fc4-8cc3-11eb-1f10-37f6d7a9d6ea
norm.(eq_velocities) # This is correct according to the guide (nm/ps)

# ╔═╡ e8fd4ae6-8cc5-11eb-3a71-6d68f6bb1191
eq_temperatures # These are way too small -- are the units incorrect?

# ╔═╡ 23f20c58-81df-11eb-010c-2bec3a114f09
begin
	Plots.plot(
		title="Temperature During Equilibration Stage",
		xlab="Time [ps]",
		ylab="Temperature [K]",
	)
	Plots.plot!(
		eq_time_range,
		eq_temperatures,
		label="simulation temperature"
	)
	Plots.plot!(
		eq_time_range,
		[reference_temp for _ ∈ 1:length(eq_time_range)],
		linestyle=:dash,
		label="reference temperature"
	)
end

# ╔═╡ 4b54bf3e-81df-11eb-17ae-b3d51c2a80d6
begin
	Plots.plot(
		eq_time_range,
		eq_energies,
		title="Energy During Equilibration Stage",
		xlab="Time [ps]",
		ylab="Energy [eV]",
		label="Total Energy"
	)
	Plots.plot!(
		eq_time_range,
		eq_kinetic_energies,
		label="Kinetic Energy"
	)
	Plots.plot!(
		eq_time_range,
		eq_potential_energies,
		label="Potential Energy"
	)
end

# ╔═╡ 3eee0850-8cd3-11eb-2519-d78c79b01022
eq_potential_energies

# ╔═╡ 7d6936a4-8cc9-11eb-000a-73547cfe8e9d
md"""
## Production Stage
"""

# ╔═╡ 075b19a8-8ccb-11eb-171a-3d8e8586716c
md"""
### Initialization
"""

# ╔═╡ 2d16a780-8cca-11eb-1768-27072238f6f0
begin
	prod_thermostat = NoThermostat()
	prod_loggers = Dict(
		"temp" => TemperatureLogger(stride),
		"coords" => CoordinateLogger(stride),
		"energy" => EnergyLogger(stride),
		"kinetic" => KineticEnergyLogger(stride),
		"potential" => PotentialEnergyLogger(stride),
		"velocities" => VelocityLogger(stride)
	)
	prod_simulation = Simulation(
		simulator=simulator,
		atoms=atoms,
		general_inters=general_inters,
		coords=deepcopy(eq_coords),
		velocities=deepcopy(eq_velocities),
		temperature=reference_temp,
		box_size=box_size,
		neighbour_finder=neighbour_finder,
		thermostat=prod_thermostat,
		loggers=prod_loggers,
		timestep=timestep,
		n_steps=prod_steps
	);
end;

# ╔═╡ 13f19d36-8ccb-11eb-0897-2db3774409b6
md"""
### Production Simulation
"""

# ╔═╡ 1ab1cd76-8ccb-11eb-3b53-4dc15ffcbb9a
begin
	simulate!(prod_simulation, parallel=true)
	prod_time_range = collect(1:prod_steps/stride) .* timestep
	prod_temperatures = prod_simulation.loggers["temp"].temperatures
	prod_energies = prod_simulation.loggers["energy"].energies
	prod_potential_energies = prod_simulation.loggers["potential"].energies
	prod_kinetic_energies = prod_simulation.loggers["kinetic"].energies
	prod_coords = prod_simulation.loggers["coords"].coords[end]
	prod_velocities = prod_simulation.loggers["velocities"].velocities[end]
end;

# ╔═╡ a716dce8-8ccb-11eb-3b9a-b38fdd95f70e
begin
	Plots.plot(
		prod_time_range,
		prod_temperatures,
		title="Temperature During Production Stage",
		xlab="Time [ps]",
		ylab="Temperature [K]",
		label="Simulation Temperature"
	)
	Plots.plot!(
		prod_time_range,
		[reference_temp for _ ∈ 1:length(prod_time_range)],
		linestyle=:dash,
		label="Reference Temperature"
	)
end

# ╔═╡ 682855b2-8ccc-11eb-1bc0-418191f1c8c6
begin
	Plots.plot(
		title="Energy During Production Stage",
		xlab="Time [ps]",
		ylab="Energy [eV]"
	)
	Plots.plot!(
		prod_time_range,
		prod_energies,
		label="Total Energy"
	)
	Plots.plot!(
		prod_time_range,
		prod_kinetic_energies,
		label="Kinetic Energy"
	)
	Plots.plot!(
		prod_time_range,
		prod_potential_energies,
		label="Potential Energy"
	)
end

# ╔═╡ d06380dc-8cc9-11eb-2181-e5233e63a8dc
md"""
## Results
"""

# ╔═╡ 2de21090-8b6d-11eb-272c-1317326981f5
begin
	initial_bins, intial_densities = rdf(initial_coords / σ, box_size)
	eq_bins, eq_densities = rdf(eq_coords / σ, box_size)
	prod_bins, prod_densities = rdf(prod_coords / σ, box_size)
	Plots.plot(	
		title="Radial Distribution Function",
		xlab="Distance r/σ",
		ylab="Radial Distribution g(r)"
	)
	Plots.plot!(
		initial_bins,
		intial_densities,
		label="initial distribution"
	)
	Plots.plot!(
		eq_bins,
		eq_densities,
		label="distribution after equilibration stage"
	)
	Plots.plot!(
		prod_bins,
		prod_densities,
		label="distribution after production stage"
	)
	# The way these are being normalized seems incorrect
end

# ╔═╡ Cell order:
# ╟─b148ddb6-817b-11eb-2878-c1b583bb35f2
# ╟─3328ea8e-8cbb-11eb-2c3f-b143741525be
# ╠═4ee38c78-8b4b-11eb-3284-45e3dd145d1f
# ╟─d8daea28-817a-11eb-26e1-a909bc31bf83
# ╟─c7559cde-817b-11eb-1d86-b5df501583b5
# ╟─baaabfde-817c-11eb-2c74-bbcff5f1d162
# ╟─cd5a9930-817c-11eb-0e02-cd5a7462b687
# ╠═d710689c-817c-11eb-0727-3b104fbbb669
# ╟─04a5c10c-8cca-11eb-1386-b128e8ce0ed8
# ╟─b4379f32-817e-11eb-2689-9d11da52ae75
# ╟─cfc54920-817e-11eb-2758-4748397ddf98
# ╠═920c3baa-817e-11eb-2037-410965f2894d
# ╟─1c301082-8180-11eb-158f-7538d52db5c3
# ╠═d2fa5508-817f-11eb-1f96-19b27221e7be
# ╠═535b18c4-8cc3-11eb-3c09-c5118d0a2d2f
# ╠═5a247fc4-8cc3-11eb-1f10-37f6d7a9d6ea
# ╠═e8fd4ae6-8cc5-11eb-3a71-6d68f6bb1191
# ╟─23f20c58-81df-11eb-010c-2bec3a114f09
# ╟─4b54bf3e-81df-11eb-17ae-b3d51c2a80d6
# ╠═3eee0850-8cd3-11eb-2519-d78c79b01022
# ╟─7d6936a4-8cc9-11eb-000a-73547cfe8e9d
# ╟─075b19a8-8ccb-11eb-171a-3d8e8586716c
# ╠═2d16a780-8cca-11eb-1768-27072238f6f0
# ╟─13f19d36-8ccb-11eb-0897-2db3774409b6
# ╠═1ab1cd76-8ccb-11eb-3b53-4dc15ffcbb9a
# ╟─a716dce8-8ccb-11eb-3b9a-b38fdd95f70e
# ╟─682855b2-8ccc-11eb-1bc0-418191f1c8c6
# ╟─d06380dc-8cc9-11eb-2181-e5233e63a8dc
# ╠═2de21090-8b6d-11eb-272c-1317326981f5
