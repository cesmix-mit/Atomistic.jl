### A Pluto.jl notebook ###
# v0.14.1

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
			"Unitful"
			"UnitfulAtomic"
			"UnitfulRecipes"
	])
	using LinearAlgebra
	using NBodySimulator
	using Plots
	using StaticArrays
	using Unitful
	using UnitfulAtomic
	using UnitfulRecipes
end

# ╔═╡ cc550040-babd-4eec-bcf0-19a687ae459d
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

# ╔═╡ a7e0d1eb-2595-4cb8-ba82-dae054c90305
begin
	M = ingredients("../src/nbs_extensions.jl");
	import .M:
		get_final_bodies
end;

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

# ╔═╡ d6697f13-a829-4126-8b3d-21f841293c1a
begin
	N = 864
	m = auconvert(6.6335209e-26u"kg")
	
	box_size = auconvert(3.47786u"nm")
	
	σ = auconvert(0.34u"nm")
	cutoff = auconvert(0.765u"nm")
	ϵ = auconvert(1.657e-21u"J")
	
	reference_temp = auconvert(94.4u"K")
	# mean_v = auconvert(√(3 * u"k" * reference_temp / m))
	mean_v = auconvert(√(u"k" * reference_temp / m))
	thermostat_prob = 0.1 # probability for Andersen Thermostat
	
	eq_steps = 2000
	prod_steps = 5000
	Δt = auconvert(1e-2u"ps")
	
	stride = 10
end;

# ╔═╡ 080622bd-6278-4bc4-b827-c53baf47d1fc
uconvert(u"m/s", mean_v)

# ╔═╡ 98cbeff2-bfed-4af6-9b67-7964a2d2d186
uconvert(u"m/s", √(3 * u"k" * reference_temp / m))

# ╔═╡ bd0903a3-e994-42a6-9650-f24473e5e39b
typeof(box_size)

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
	initial_bodies = generate_bodies_in_cell_nodes(
		N, ustrip(m), ustrip(mean_v), ustrip(box_size)
	)
	potentials = Dict(
		:lennard_jones => LennardJonesParameters(ustrip(ϵ), ustrip(σ), ustrip(cutoff))
	)
	eq_system = PotentialNBodySystem(initial_bodies, potentials)
	
	boundary_conditions = CubicPeriodicBoundaryConditions(ustrip(box_size))
	eq_thermostat = AndersenThermostat(
		ustrip(reference_temp), thermostat_prob / ustrip(Δt)
	)
	eq_simulation = NBodySimulation(
		eq_system,
		(0.0, eq_steps * ustrip(Δt)),
		boundary_conditions,
		eq_thermostat,
		1.0
	)
	
	simulator = VelocityVerlet()
end;

# ╔═╡ 64b0e53c-91f4-11eb-0377-3fe1fd44bfef
md"""
### Equilibration Simulation
"""

# ╔═╡ 673c960c-91f4-11eb-1e26-f543ed4c7f5b
begin
	eq_result = run_simulation(eq_simulation, simulator, dt=ustrip(Δt))
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
		(eq_steps * ustrip(Δt), (eq_steps + prod_steps) * ustrip(Δt)),
		boundary_conditions,
		1.0
	)
end;

# ╔═╡ 7d3ada68-91f4-11eb-0ef0-a750182f882c
md"""
### Production Simulation
"""

# ╔═╡ 895ff2e2-91f4-11eb-2610-f3aea31b271f
begin
	prod_result = run_simulation(prod_simulation, simulator, dt=ustrip(Δt));
	prod_bodies = get_final_bodies(prod_result)
end;

# ╔═╡ 8a9178c0-91f4-11eb-36f1-61f957fe5f8f
md"""
## Results
"""

# ╔═╡ 2d571c96-9251-11eb-062f-9d0262fc3ad6
begin
	eq_time_range = [auconvert(u"ps", t) for (i,t) ∈ enumerate(eq_result.solution.t) if (i - 1) % stride == 0]
	prod_time_range = [auconvert(u"ps", t) for (i,t) ∈ enumerate(prod_result.solution.t) if (i - 1) % stride == 0]
	total_time_range = vcat(eq_time_range, prod_time_range)
end;

# ╔═╡ 89fc84ca-d6ec-4847-8c04-a4028494f28c
auconvert(u"K", temperature(eq_result, 0))

# ╔═╡ 91061d14-91f4-11eb-306c-83c482fdf378
begin
	plot(
		title="Temperature during Simulation",
		xlab="Time",
		ylab="Temperature",
	)
	plot!(
		eq_time_range,
		t -> auconvert(u"K", temperature(eq_result, ustrip(auconvert(t)))),
		label="Equilibrium Stage Temperature",
		color=1
	)
	plot!(
		prod_time_range,
		t -> auconvert(u"K", temperature(prod_result, ustrip(auconvert(t)))),
		label="Production Stage Temperature",
		color=3
	)
	plot!(
		total_time_range,
		t -> uconvert(u"K", reference_temp),
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
		xlab="Time",
		ylab="Energy",
		legend=:right
	)
	plot!(
		eq_time_range,
		t -> auconvert(u"hartree", kinetic_energy(eq_result, ustrip(auconvert(t)))),
		label="Kinetic Energy",
		color=2
	)
	plot!(
		eq_time_range,
		t -> auconvert(u"hartree", potential_energy(eq_result, ustrip(auconvert(t)))),
		label="Potential Energy",
		color=1
	)
	plot!(
		eq_time_range,
		t -> auconvert(u"hartree", total_energy(eq_result, ustrip(auconvert(t)))),
		label="Total Energy",
		color=3
	)
	plot!(
		prod_time_range,
		t -> auconvert(u"hartree", kinetic_energy(prod_result, ustrip(auconvert(t)))),
		label=false,
		color=2
	)
	plot!(
		prod_time_range,
		t -> auconvert(u"hartree", potential_energy(prod_result, ustrip(auconvert(t)))),
		label=false,
		color=1
	)
	plot!(
		prod_time_range,
		t -> auconvert(u"hartree", total_energy(prod_result, ustrip(auconvert(t)))),
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
# end;

# ╔═╡ 656a0a6c-9253-11eb-23f4-4bddb7746a4d
# begin
# 	plot(
# 		title="Radial Distribution Function",
# 		xlab="Distance r/σ",
# 		ylab="Radial Distribution"
# 	)
# 	plot!(
# 		[auconvert(u"bohr", r) for r ∈ eq_rs] / σ,
# 		eq_grf,
# 		label="Equilibrium Stage Distribution"
# 	)
# 	plot!(
# 		[auconvert(u"bohr", r) for r ∈ prod_rs] / σ,
# 		prod_grf,
# 		label="Production Stage Distribution"
# 	)
# end

# ╔═╡ Cell order:
# ╟─d0806d6e-91f2-11eb-045b-fd5b1b63acb1
# ╟─cc550040-babd-4eec-bcf0-19a687ae459d
# ╟─a7e0d1eb-2595-4cb8-ba82-dae054c90305
# ╟─35f8a14a-91f3-11eb-2cac-6da99c7f96ac
# ╟─1b675294-91f4-11eb-1ad7-a5ca38362fb2
# ╟─28bf50ae-91f4-11eb-1581-dba8ca5e5b76
# ╟─3464c33c-91f4-11eb-2507-9db654e50fb3
# ╠═d6697f13-a829-4126-8b3d-21f841293c1a
# ╠═080622bd-6278-4bc4-b827-c53baf47d1fc
# ╠═98cbeff2-bfed-4af6-9b67-7964a2d2d186
# ╠═bd0903a3-e994-42a6-9650-f24473e5e39b
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
# ╠═89fc84ca-d6ec-4847-8c04-a4028494f28c
# ╟─91061d14-91f4-11eb-306c-83c482fdf378
# ╟─ec2ba49e-9243-11eb-02ab-0dea8982f775
# ╠═8646dad8-9246-11eb-36c8-f32871baa335
# ╠═656a0a6c-9253-11eb-23f4-4bddb7746a4d
