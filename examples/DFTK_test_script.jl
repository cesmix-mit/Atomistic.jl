using Atomistic
using DFTK
using Unitful
using UnitfulAtomic

box_size = 1.5u"nm"
lattice = box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]]

psp = ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier))
atoms = [psp => [
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.75],
    [0.25, 0.75, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.25],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
    [0.75, 0.75, 0.75],
]]

model = model_LDA(lattice, atoms)
Ecut = 5u"hartree"
kgrid = [1, 1, 1]

PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)