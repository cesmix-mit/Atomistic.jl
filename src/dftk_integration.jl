using DFTK
using NBodySimulator
using StaticArrays
using Unitful
using UnitfulAtomic

include("./molecular_simulation.jl")

setup_threading()

struct DFTKPotentialParameters{pType <: Real} <: PotentialParameters
	forces::Vector{SVector{3, pType}}
end

function NBodySimulator.get_accelerating_function(parameters::DFTKPotentialParameters, simulation::NBodySimulation)
    masses = get_masses(simulation.system)
    (dv, u, v, t, i) -> begin dv .+= parameters.forces[i] / masses[i] end
end

struct DFTKForceGenerationParameters <: ForceGenerationParameters
    box_size::Quantity
    lattice::AbstractArray{Quantity, 2}
    psp::ElementPsp
    kgrid::AbstractVector{Integer}
    Ecut::Quantity
    tolerance::AbstractFloat
end

function calculate_scf(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters)
    atoms = [parameters.psp => [auconvert.(u"bohr", b.r) / parameters.box_size for b ∈ bodies]]

    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    return @time self_consistent_field(basis, tol=parameters.tolerance)
end

function generate_forces(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters)
    scfres = calculate_scf(bodies, parameters)
    return DFTKPotentialParameters(compute_forces_cart(scfres)[1])
end

function analyze_convergence(bodies::Vector{MassBody}, parameters::DFTKForceGenerationParameters, cutoffs::Vector{<:Quantity})
    options = Dict(Ecut => DFTKForceGenerationParameters(parameters.box_size, parameters.lattice, parameters.psp, parameters.kgrid, Ecut, parameters.tolerance) for Ecut in cutoffs)
    fields = Dict(Ecut => calculate_scf(bodies, options[Ecut]) for Ecut in cutoffs)
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