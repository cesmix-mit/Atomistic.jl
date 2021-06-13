using DFTK
using LinearAlgebra

setup_threading()

lattice = 20diagm([1,1,1.])
Ar = ElementPsp(:Ar, psp=load_psp("hgh/lda/ar-q8.hgh"))
pos = map(v -> lattice \ v + 0.5ones(3), vec([[3i, 3j, 3k] for i in 0:1, j in 0:1, k in 0:1]))
atoms = [Ar => pos]
Ecut = 10
model = model_LDA(lattice, atoms, temperature=1e-3)
basis = PlaneWaveBasis(model, Ecut, kgrid=[1, 1, 1])

DFTK.reset_timer!(DFTK.timer)
scfres = self_consistent_field(basis, tol=1e-6, α=0.7, mixing=LdosMixing())

println(scfres.energies)
println()
println()
forces = compute_forces_cart(scfres)
println(only(forces))
println()
println()
println(DFTK.timer)