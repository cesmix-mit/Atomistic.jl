using DFTK
using Plots
using Unitful
using UnitfulAtomic

### Copied from https://juliamolsim.github.io/DFTK.jl/stable/guide/tutorial/

println("step 1")
# 1. Define lattice and atomic positions
a = 5.431u"angstrom"          # Silicon lattice constant
lattice = a / 2 * [[0 1 1.];  # Silicon lattice vectors
                   [1 0 1.];  # specified column by column
                   [1 1 0.]]

# Load HGH pseudopotential for Silicon
Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))

# Specify type and positions of atoms
atoms = [Si => [ones(3) / 8, -ones(3) / 8]]

println("step 2")
# 2. Select model and basis
model = model_LDA(lattice, atoms)
kgrid = [4, 4, 4]     # k-point grid (Regular Monkhorst-Pack grid)
Ecut = 7              # kinetic energy cutoff
# Ecut = 190.5u"eV"  # Could also use eV or other energy-compatible units
basis = PlaneWaveBasis(model, Ecut; kgrid=kgrid)

println("step 3")
# 3. Run the SCF procedure to obtain the ground state
@time scfres = self_consistent_field(basis, tol=1e-8);