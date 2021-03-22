### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 2cc0ec2e-7b19-11eb-082e-5d379b88f087
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			"Makie"
			"Molly"
			"Plots"
			"PlutoUI"
	])
	using Molly
	using Makie
	using Plots
	using PlutoUI
end

# ╔═╡ b7ce7172-7b19-11eb-1c9b-49731d0c8da3
md"""
# Exploration of Molly.jl
Documentation can be referenced [here](https://juliamolsim.github.io/Molly.jl/stable/).
"""

# ╔═╡ 6ef69b08-7c3d-11eb-195b-479b53f9d132
md"""
## Initial Example
"""

# ╔═╡ 9a934f2c-7c3d-11eb-3ef5-ad296067cb16
md"""
Below is the example given in the documentation for simulating a Lennard-Jones gas which I will break into steps to explain.
"""

# ╔═╡ fd53322a-7c3e-11eb-1da9-5963abc0c560
md"""
First we must define a list of atoms that we are working with. Each atom is defined by 
- an electric charge `charge`
- a mass `mass`
- the Lennard-Jones finite distance at which the inter-particle potential is zero `σ`
- the Lennard-Jones depth of the potential well `ϵ`
In this example, we are setting all atoms to be uncharged.
"""

# ╔═╡ daaf500a-7c3e-11eb-1538-f1fcddb9c12c
begin
	n_atoms = 100
	mass = 10.0
	atoms = [Atom(mass=mass, σ=0.3, ϵ=0.2) for i in 1:n_atoms]
end

# ╔═╡ b0f83a00-7c3f-11eb-1e0d-cbcf46dd3192
md"""
Then we also need initial conditions for each atom's position and velocity. Static vectors in 3-space are used for performance. The `velocity` function generates a random velocity from the Maxwell-Boltzmann distribution.
"""

# ╔═╡ c6b9938e-7c3f-11eb-3fdc-2fecb3c54175
begin
	box_size = 2.0 # nm
	coords = [box_size .* rand(SVector{3}) for i in 1:n_atoms]
end

# ╔═╡ 6870f2a8-7c45-11eb-1a3b-bd2f7ddddbe4
begin
	temp = 100 # K
	velocities = [velocity(mass, temp) for i in 1:n_atoms]
end

# ╔═╡ 2142353a-7c41-11eb-25ad-9de746eb598a
md"""
We also define a tuple of interactions between the atoms. In this case we are using the built-in [Lennard-Jones potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) which is defined by the function

$V_{LJ} = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-\left(\frac{\sigma}{r}\right)^{6}\right]$
"""

# ╔═╡ 46c99e24-7c41-11eb-2779-67978c7ed4df
general_inters = (LennardJones(),)

# ╔═╡ 8321eeca-7c42-11eb-2edd-897caa16dc7c
md"""
The measurements during the simulation will be determined by the loggers that we use. In this case we are measuring the temperature and coordinates after every 10 timesteps.
"""

# ╔═╡ e711615e-7c42-11eb-1a25-0b0da724b893
begin
	sample_rate = 10
	loggers = Dict(
		"temp" => TemperatureLogger(sample_rate),
		"coords" => CoordinateLogger(sample_rate)
	)
end

# ╔═╡ 2e344ac0-7c42-11eb-0865-f1e3068250d2
md"""
We can now create the simulation object. In this case we are simulating over 1,000 0.002 ps timesteps using the [velocity Verlet algorithm](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet). We also use an [Andersen thermostat](https://en.wikipedia.org/wiki/Andersen_thermostat) to keep a constant termperature.
"""

# ╔═╡ 879c5396-7c42-11eb-2743-19c295d57281
begin
	timestep = 0.002 # ps
    n_steps = 1_000
	
	s = Simulation(
		simulator=VelocityVerlet(), # Use velocity Verlet integration
		atoms=atoms,
		general_inters=general_inters,
		coords=coords,
		velocities=velocities,
		temperature=temp,
		box_size=box_size,
		thermostat=AndersenThermostat(1.0), # Coupling constant of 1.0
		loggers=loggers,
		timestep=timestep,
		n_steps=n_steps
	)
end

# ╔═╡ b4168990-7c43-11eb-2e27-95a5f656aeb5
md"""
The `simulate!` method runs the simulation. We can use the `parallel` argument to run the simulation on all available threads. This method returns the coordinates of the atoms at the end of the simulation.
"""

# ╔═╡ 71ecb86e-7c43-11eb-0453-1d351b2bb5eb
simulate!(s, parallel=true)

# ╔═╡ f963d2aa-7c43-11eb-0d32-95f9338361c4
s.loggers["temp"].temperatures

# ╔═╡ 1bbee6c8-7c5a-11eb-257c-9b952a96d207
Plots.plot(s.loggers["temp"].temperatures, title="Temperature During Simulation", xlab="Sample Number ($(sample_rate * timestep)ps Increments)", ylab="Temperature [K]", legend=false)

# ╔═╡ 3a631c3a-7c45-11eb-0801-91f5abc6bfac
s.loggers["coords"].coords

# ╔═╡ 31dedb40-7c47-11eb-3bc4-3949dffaf1d9
md"""
We can use [Makie.jl](http://makie.juliaplots.org/stable/) to visualize the simulation with the `visualize` function.
"""

# ╔═╡ 4697498c-7c47-11eb-11e0-dd4c1c40b300
file = visualize(s.loggers["coords"], box_size, "../artifacts/example_simulation.gif", markersize=10);

# ╔═╡ 2ae89586-7c53-11eb-05fd-1377defa0aab
md"""$(LocalResource(file))"""

# ╔═╡ 781a77ee-7c4d-11eb-16ec-31a1e4773418
md"""
## Investigating Potentials
"""

# ╔═╡ bbfb64fe-7c53-11eb-27a9-7dc6e2e7e444
md"""
The forces between particles are derived from the underlying potential function that is being used.

$$\vec{F_i}=-\sum_j\frac{dV_{ij}(r_{ij})}{dr_{ij}}\frac{\vec{r_{ij}}}{r_{ij}}$$

Molly defines `GeneralInteractions` that act between all or most atoms and `SpecificInteractions` that act between specific atoms.

The built-in general interactions are `LennardJones`, `SoftSphere`, `Mie`, `Coulomb`, and `Gravity`.

Custom general interactions can be definined by subtyping `GeneralInteraction`. This involves
- setting a boolean value `nl_only` which determines whether only neighboring atoms or all atoms in the simulation are considered
- defining a force method that governs how two particles interact

The built-in specific interactions are `HarmonicBond`, `HarmonicAngle`, and `Torsion`.

Custom specific interactions can be defined by subtyping `SpecificInteraction`. This involves
- storing references to the particular particles that exhibit the interaction
- defining a force methods that governs how these specific particles interact
"""

# ╔═╡ 69ff4a44-7c53-11eb-048d-cd7fd21ad475
LennardJones()

# ╔═╡ fac4f1bc-7c55-11eb-0483-53887c03b1f4
md"""
## Other Configuration Options and Notes
"""

# ╔═╡ 0b7c4642-7c57-11eb-1afa-5bacac28c735
md"""
Molly supports the `VelocityVerlet` and `VelocityFreeVerlet` out of the box and we can implement other simulators by subtyping `Simulator`.

The built-in loggers are `TemperatureLogger`, `CoordinateLogger`, `EnergyLogger`, and `StructureWriter` and we can also implemenet additional loggers by subtyping `Logger`.

There is also support for customizing the `Thermostat` and `NeighborhoodFinder`.
"""

# ╔═╡ 7ae0d3ba-7c57-11eb-35b0-150d86d14c1a
md"""
There are examples in the documentation of importing `.top` and `.gro` files in order to initialize the simulation. This is done through the `readinputs` method which takes in a `topology_file` and a `coordinate_file`. I am not familiar with these file formats, but they might be useful.
"""

# ╔═╡ be08600a-7c56-11eb-04d0-b31354a4e783
md"""
## GPU Acceleration
"""

# ╔═╡ 17786cc2-7c57-11eb-351b-055682a1dc54
md"""
Simulations can be run on the GPU of a CUDA-compatible device. Doing this requires replacing the `Atom` type with the `AtomMin` type, and the documentation also recommends using `Float32`s when running on the GPU. I haven't tested this out yet, but it's something we can look into more in the future.
"""

# ╔═╡ Cell order:
# ╟─2cc0ec2e-7b19-11eb-082e-5d379b88f087
# ╟─b7ce7172-7b19-11eb-1c9b-49731d0c8da3
# ╟─6ef69b08-7c3d-11eb-195b-479b53f9d132
# ╟─9a934f2c-7c3d-11eb-3ef5-ad296067cb16
# ╟─fd53322a-7c3e-11eb-1da9-5963abc0c560
# ╠═daaf500a-7c3e-11eb-1538-f1fcddb9c12c
# ╟─b0f83a00-7c3f-11eb-1e0d-cbcf46dd3192
# ╠═c6b9938e-7c3f-11eb-3fdc-2fecb3c54175
# ╠═6870f2a8-7c45-11eb-1a3b-bd2f7ddddbe4
# ╟─2142353a-7c41-11eb-25ad-9de746eb598a
# ╠═46c99e24-7c41-11eb-2779-67978c7ed4df
# ╟─8321eeca-7c42-11eb-2edd-897caa16dc7c
# ╠═e711615e-7c42-11eb-1a25-0b0da724b893
# ╟─2e344ac0-7c42-11eb-0865-f1e3068250d2
# ╠═879c5396-7c42-11eb-2743-19c295d57281
# ╟─b4168990-7c43-11eb-2e27-95a5f656aeb5
# ╠═71ecb86e-7c43-11eb-0453-1d351b2bb5eb
# ╠═f963d2aa-7c43-11eb-0d32-95f9338361c4
# ╠═1bbee6c8-7c5a-11eb-257c-9b952a96d207
# ╠═3a631c3a-7c45-11eb-0801-91f5abc6bfac
# ╟─31dedb40-7c47-11eb-3bc4-3949dffaf1d9
# ╠═4697498c-7c47-11eb-11e0-dd4c1c40b300
# ╟─2ae89586-7c53-11eb-05fd-1377defa0aab
# ╟─781a77ee-7c4d-11eb-16ec-31a1e4773418
# ╟─bbfb64fe-7c53-11eb-27a9-7dc6e2e7e444
# ╠═69ff4a44-7c53-11eb-048d-cd7fd21ad475
# ╟─fac4f1bc-7c55-11eb-0483-53887c03b1f4
# ╟─0b7c4642-7c57-11eb-1afa-5bacac28c735
# ╟─7ae0d3ba-7c57-11eb-35b0-150d86d14c1a
# ╟─be08600a-7c56-11eb-04d0-b31354a4e783
# ╟─17786cc2-7c57-11eb-351b-055682a1dc54
