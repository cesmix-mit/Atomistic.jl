using DFTK
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

include("./molecular_simulation.jl")

setup_threading()

struct DFTKForceGenerationParameters <: ForceGenerationParameters
    lattice::AbstractArray{Quantity, 2}
    psp::ElementPsp
    kgrid::AbstractVector{Integer}
    Ecut::Quantity
    tolerance::AbstractFloat
end

struct DFTKPotentialParameters{pType <: Real} <: PotentialParameters
	forces::Vector{SVector{3, pType}}
end

# import NBodySimulator.get_accelerating_function
function get_accelerating_function(parameters::DFTKPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin dv .+= parameters.forces[i] / masses[i] end
end

function calculate_scf(bodies::Vector{MassBody}, box_size::Quantity, parameters::DFTKForceGenerationParameters)
    atoms = [parameters.psp => [auconvert.(u"bohr", b.r) / box_size for b ∈ bodies]]

    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    @time scfres = self_consistent_field(basis, tol=parameters.tolerance)
end

function generate_forces(bodies::Vector{MassBody}, box_size::Quantity, parameters::DFTKForceGenerationParameters)
    scfres = calculate_scf(bodies, box_size, parameters)
    return DFTKPotentialParameters(compute_forces_cart(scfres)[1])
end

function analyze_convergence(bodies::Vector{MassBody}, box_size::Quantity, parameters::DFTKForceGenerationParameters, cutoffs::Vector{<:Quantity})
    options = Dict(Ecut => DFTKForceGenerationParameters(parameters.lattice, parameters.psp, parameters.kgrid, Ecut, parameters.tolerance) for Ecut in cutoffs)
    fields = Dict(Ecut => calculate_scf(bodies, box_size, options[Ecut]) for Ecut in cutoffs)
    energies = Dict(Ecut => fields[Ecut].energies.total for Ecut in cutoffs)
    plot(
		title="DFTK Analysis",
		xlab="Ecut",
		ylab="Total Energy",
        legend=false
	)
	plot!(
        cutoffs,
        c -> auconvert(u"hartree", energies[c])
	)
end