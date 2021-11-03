using Atomistic    
using DFTK    
using Unitful    
using UnitfulAtomic    
disable_threading()    
    
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
    
model = model_LDA(lattice, atoms; symmetries=false)    
Ecut = 5u"hartree"    
kgrid = [2, 1, 1]    
    
basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)    
mpi_master() && display(basis)    
    
scfres = self_consistent_field(basis)    
    
mpi_master() && display(scfres.energies) 