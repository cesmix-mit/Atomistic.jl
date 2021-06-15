setup_threading(n_blas=4)
using PyCall

Base.@kwdef struct DFTKForceGenerationParameters <: ForceGenerationParameters
    box_size::Quantity
    psp::ElementPsp
    lattice::AbstractArray{Quantity, 2}
    Ecut::Quantity
    kgrid::AbstractVector{Integer}
    n_bands::Union{Integer, Missing} = missing
    tol::Union{AbstractFloat, Missing} = missing
    α::Union{AbstractFloat, Missing} = missing
    mixing = missing # There is no abstract type for mixing :(
    previousscfres = Any[nothing]
end
DFTKForceGenerationParameters(parameters::DFTKForceGenerationParameters, Ecut::Quantity) = DFTKForceGenerationParameters(parameters.box_size, parameters.psp, parameters.lattice, Ecut, parameters.kgrid, parameters.n_bands, parameters.tol, parameters.α, parameters.mixing);

function calculate_scf(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters)
    function trunc(r)
        r - floor(Int, r)
    end
    atoms = [parameters.psp => [trunc.(b.r / austrip(parameters.box_size)) for b ∈ bodies]]

    asa = ase_atoms(austrip.(parameters.lattice), atoms)
    traj = pyimport("ase.io.trajectory").Trajectory("run.traj", "a", asa)
    traj.write(asa)
    traj.close()

    model = model_LDA(parameters.lattice, atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    extraargs = (; )
    if parameters.previousscfres[1] !== nothing
        scfres = parameters.previousscfres[1]
        extraargs = (ψ=scfres.ψ, ρ=scfres.ρ, )
    end

    scfres =  @time self_consistent_field(basis; extraargs..., (f=>getfield(parameters, f) for f in (:n_bands, :tol, :α, :mixing) if getfield(parameters, f) !== missing)...)

    parameters.previousscfres[1] = scfres
    scfres
end

function generate_forces(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters)
    scfres = calculate_scf(bodies, parameters)
    return compute_forces_cart(scfres)[1]
end

function analyze_convergence(bodies::AbstractVector{MassBody}, parameters::DFTKForceGenerationParameters, cutoffs::AbstractVector{<:Quantity})
    options = Dict(Ecut => DFTKForceGenerationParameters(parameters, Ecut) for Ecut in cutoffs)
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
