using DFTK
using Unitful
using UnitfulAtomic

setup_threading()

function calculate_forces(bodies::Array{MassBody}, box_size::Quantity)
    # 1. Define lattice and atomic positions
    # a = 5.26u"angstrom"          # Argon lattice constant
    lattice = box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]]

    # Load HGH pseudopotential for Argon
    Ar = ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier))

    # Specify type and positions of atoms
    atoms = [Ar => [auconvert.(u"bohr", b.r) / box_size for b ∈ bodies]]

    # 2. Select model and basis
    model = model_LDA(lattice, atoms)
    kgrid = [1, 1, 1]     # k-point grid (Regular Monkhorst-Pack grid)
    Ecut = 10             # kinetic energy cutoff
    # Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
    basis = PlaneWaveBasis(model, Ecut; kgrid=kgrid)

    # 3. Run the SCF procedure to obtain the ground state
    @time scfres = self_consistent_field(basis, tol=1e-4)

    return compute_forces_cart(scfres)[1]
end