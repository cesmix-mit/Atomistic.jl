### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ d0806d6e-91f2-11eb-045b-fd5b1b63acb1
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			"LinearAlgebra"
			"NBodySimulator"
			"Plots"
			"StaticArrays"
	])
	using LinearAlgebra
	using NBodySimulator
	using Plots
	using StaticArrays
end

# ╔═╡ 35f8a14a-91f3-11eb-2cac-6da99c7f96ac
md"""
# Example MD Simulation of Liquid Argon using NBodySimulator
"""

# ╔═╡ 1b675294-91f4-11eb-1ad7-a5ca38362fb2
md"""
This is a basic simulation of a cluster of Argon molecules in a cubical box. The purpose of this notebook is to experiment with NBodySimulator.jl with the built-in Lennard-Jones potential and demonstrate how to use the package to reproduce a basic known experiment.
"""

# ╔═╡ 28bf50ae-91f4-11eb-1581-dba8ca5e5b76
md"""
## Parameters
"""

# ╔═╡ 3464c33c-91f4-11eb-2507-9db654e50fb3
md"""
The parameters for this example are from [this guide](https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf).
"""

# ╔═╡ 5827e012-9245-11eb-0eff-39184736df95
begin
	J_per_eV = 1.602176634e-19
	kg_per_amu = 1.66054e-27
end;

# ╔═╡ 44b6546a-91f4-11eb-15bd-514d01ef2b52
begin
	kb = 1.380649e-23 # J/K
	
	N = 864
	m = 6.6335209e-26 # kg
	
	box_size = 3.47786e-9 # m
	
	σ = 0.34e-9 # m
	cutoff = 0.765e-9 # m
	ϵ = 1.657e-21 # J
	
	reference_temp = 94.4 # K
	mean_v = √(kb * reference_temp / m) # m/s
	thermostat_prob = 0.01
	
	eq_steps = 2000
	prod_steps = 5000
	Δt = 1e-14 # s
	
	stride = 10
end;

# ╔═╡ e5609e6a-924a-11eb-22c0-d14c09e79523
md"""
## Utility Functions
"""

# ╔═╡ f0b571b4-924a-11eb-0f3b-95b8fe895e5c
function get_final_bodies(sr::NBodySimulator.SimulationResult)
	positions = get_position(sr, sr.solution.t[end])
	velocities = get_velocity(sr, sr.solution.t[end])
	masses = get_masses(sr.simulation.system)
	N = length(masses)
	bodies = MassBody[]
	for i ∈ 1:N
		push!(bodies, MassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i]))
	end
	return bodies
end

# ╔═╡ 4b5184ac-91f4-11eb-1dd7-b537e60cbe69
md"""
## Equilibration Stage
"""

# ╔═╡ 539a6624-91f4-11eb-2008-ad73091ea78d
md"""
The referenced guide uses a GROMACS coordinate file for the initialization, but we are just using a random distribution over the volume of the box.
"""

# ╔═╡ 097cf250-9241-11eb-3489-c70ca430c691
md"""
### Initialization
"""

# ╔═╡ 5aa365e4-91f4-11eb-291e-79be718a7132
begin
	initial_bodies = generate_bodies_in_cell_nodes(N, m, mean_v, box_size)
	potentials = Dict(:lennard_jones => LennardJonesParameters(ϵ, σ, cutoff))
	eq_system = PotentialNBodySystem(initial_bodies, potentials)
	
	boundary_conditions = CubicPeriodicBoundaryConditions(box_size)
	eq_thermostat = AndersenThermostat(reference_temp, thermostat_prob / Δt)
	eq_simulation = NBodySimulation(
		eq_system,
		(0.0, eq_steps * Δt),
		boundary_conditions,
		eq_thermostat,
		kb
	)
	
	simulator = VelocityVerlet()
end;

# ╔═╡ 64b0e53c-91f4-11eb-0377-3fe1fd44bfef
md"""
### Equilibration Simulation
"""

# ╔═╡ 673c960c-91f4-11eb-1e26-f543ed4c7f5b
begin
	eq_result = run_simulation(eq_simulation, simulator, dt=Δt)
	eq_bodies = get_final_bodies(eq_result)
end;

# ╔═╡ 6ec3d124-91f4-11eb-0b1b-01a3a56c160e
md"""
## Production Stage
"""

# ╔═╡ 70396e06-91f4-11eb-26c0-3940a8f546af
md"""
### Initialization
"""

# ╔═╡ 71a390b4-91f4-11eb-2a16-492dffa48054
begin
	prod_system = PotentialNBodySystem(eq_bodies, potentials)
	
	prod_simulation = NBodySimulation(
		prod_system,
		(eq_steps * Δt, (eq_steps + prod_steps) * Δt),
		boundary_conditions,
		kb
	)
end;

# ╔═╡ 7d3ada68-91f4-11eb-0ef0-a750182f882c
md"""
### Production Simulation
"""

# ╔═╡ 895ff2e2-91f4-11eb-2610-f3aea31b271f
begin
	prod_result = run_simulation(prod_simulation, simulator, dt=Δt);
	prod_bodies = get_final_bodies(prod_result)
end;

# ╔═╡ 8a9178c0-91f4-11eb-36f1-61f957fe5f8f
md"""
## Results
"""

# ╔═╡ 2d571c96-9251-11eb-062f-9d0262fc3ad6
begin
	eq_time_range = [t * 1e12 for (i,t) ∈ enumerate(eq_result.solution.t) if (i - 1) % stride == 0]
	prod_time_range = [t * 1e12 for (i,t) ∈ enumerate(prod_result.solution.t) if (i - 1) % stride == 0]
	total_time_range = vcat(eq_time_range, prod_time_range)
end;

# ╔═╡ 91061d14-91f4-11eb-306c-83c482fdf378
begin
	plot(
		title="Temperature during Simulation",
		xlab="Time [ps]",
		ylab="Temperature [K]",
	)
	plot!(
		eq_time_range,
		t -> temperature(eq_result, t * 1e-12),
		label="Equilibrium Stage Temperature",
		color=1
	)
	plot!(
		prod_time_range,
		t -> temperature(prod_result, t * 1e-12),
		label="Production Stage Temperature",
		color=3
	)
	plot!(
		total_time_range,
		t -> reference_temp,
		label="Reference Temperature",
		color=2,
		linestyle=:dash
	)
	vline!(
		[eq_time_range[end]],
		label=false,
		color=:black,
		linestyle=:dash
	)
end

# ╔═╡ ec2ba49e-9243-11eb-02ab-0dea8982f775
begin
	plot(
		title="Energy during Simulation",
		xlab="Time [ps]",
		ylab="Energy [eV]",
		legend=:right
	)
	plot!(
		eq_time_range,
		t -> kinetic_energy(eq_result, t * 1e-12) / J_per_eV,
		label="Kinetic Energy",
		color=2
	)
	plot!(
		eq_time_range,
		t -> potential_energy(eq_result, t * 1e-12) / J_per_eV,
		label="Potential Energy",
		color=1
	)
	plot!(
		eq_time_range,
		t -> total_energy(eq_result, t * 1e-12) / J_per_eV,
		label="Total Energy",
		color=3
	)
	plot!(
		prod_time_range,
		t -> kinetic_energy(prod_result, t * 1e-12) / J_per_eV,
		label=false,
		color=2
	)
	plot!(
		prod_time_range,
		t -> potential_energy(prod_result, t * 1e-12) / J_per_eV,
		label=false,
		color=1
	)
	plot!(
		prod_time_range,
		t -> total_energy(prod_result, t * 1e-12) / J_per_eV,
		label=false,
		color=3
	)
	vline!(
		[eq_time_range[end]],
		label=false,
		color=:black,
		linestyle=:dash
	)
end

# ╔═╡ 8646dad8-9246-11eb-36c8-f32871baa335
# begin
# 	eq_rs, eq_grf = rdf(eq_result)
# 	prod_rs, prod_grf = rdf(prod_result)
# end

# ╔═╡ 656a0a6c-9253-11eb-23f4-4bddb7746a4d
# begin
# 	plot(
# 		title="Radial Distribution Function",
# 		xlab="Distance r/σ",
# 		ylab="Radial Distribution g(r)"
# 	)
# 	plot!(
# 		eq_rs / σ,
# 		eq_grf,
# 		label="Equilibrium Stage Distribution"
# 	)
# 	plot!(
# 		prod_rs / σ,
# 		prod_grf,
# 		label="Production Stage Distribution"
# 	)
# end

# ╔═╡ Cell order:
# ╟─d0806d6e-91f2-11eb-045b-fd5b1b63acb1
# ╟─35f8a14a-91f3-11eb-2cac-6da99c7f96ac
# ╟─1b675294-91f4-11eb-1ad7-a5ca38362fb2
# ╟─28bf50ae-91f4-11eb-1581-dba8ca5e5b76
# ╟─3464c33c-91f4-11eb-2507-9db654e50fb3
# ╠═5827e012-9245-11eb-0eff-39184736df95
# ╠═44b6546a-91f4-11eb-15bd-514d01ef2b52
# ╟─e5609e6a-924a-11eb-22c0-d14c09e79523
# ╟─f0b571b4-924a-11eb-0f3b-95b8fe895e5c
# ╟─4b5184ac-91f4-11eb-1dd7-b537e60cbe69
# ╟─539a6624-91f4-11eb-2008-ad73091ea78d
# ╟─097cf250-9241-11eb-3489-c70ca430c691
# ╠═5aa365e4-91f4-11eb-291e-79be718a7132
# ╟─64b0e53c-91f4-11eb-0377-3fe1fd44bfef
# ╠═673c960c-91f4-11eb-1e26-f543ed4c7f5b
# ╟─6ec3d124-91f4-11eb-0b1b-01a3a56c160e
# ╟─70396e06-91f4-11eb-26c0-3940a8f546af
# ╠═71a390b4-91f4-11eb-2a16-492dffa48054
# ╟─7d3ada68-91f4-11eb-0ef0-a750182f882c
# ╠═895ff2e2-91f4-11eb-2610-f3aea31b271f
# ╟─8a9178c0-91f4-11eb-36f1-61f957fe5f8f
# ╠═2d571c96-9251-11eb-062f-9d0262fc3ad6
# ╟─91061d14-91f4-11eb-306c-83c482fdf378
# ╟─ec2ba49e-9243-11eb-02ab-0dea8982f775
# ╠═8646dad8-9246-11eb-36c8-f32871baa335
# ╠═656a0a6c-9253-11eb-23f4-4bddb7746a4d
