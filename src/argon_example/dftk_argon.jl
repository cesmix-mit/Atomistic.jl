using DFTK
using Unitful
using UnitfulAtomic

setup_threading()

function argon_scf(bodies::Vector{MassBody}, box_size::Quantity, Ecut::Quantity)
    # 1. Define lattice and atomic positions
    lattice = box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]]

    # Load HGH pseudopotential for Argon
    Ar = ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier))

    # Specify type and positions of atoms
    atoms = [Ar => [auconvert.(u"bohr", b.r) / box_size for b ∈ bodies]]

    # 2. Select model and basis
    model = model_LDA(lattice, atoms)
    kgrid = [1, 1, 1]     # k-point grid (Regular Monkhorst-Pack grid)
    basis = PlaneWaveBasis(model, Ecut; kgrid=kgrid)

    # 3. Run the SCF procedure to obtain the ground state
    return @time scfres = self_consistent_field(basis, tol=1e-4)
end

function calculate_dftk_forces(bodies::Vector{MassBody}, box_size::Quantity, Ecut::Quantity=10u"hartree")
    scfres = argon_scf(bodies, box_size, Ecut)
    return [auconvert.(u"hartree/bohr", f) for f in compute_forces_cart(scfres)[1]]
end

function analyze_convergence(bodies::Vector{MassBody}, box_size::Quantity, cutoffs::Vector=[e * u"hartree" for e in (10, 12, 14, 16, 18, 20)])
    fields = Dict(Ecut => argon_scf(bodies, box_size, Ecut) for Ecut in cutoffs)
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