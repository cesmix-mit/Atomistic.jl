# TODO: A bit fragile because I'm assuming the bare minimum set of atoms_data keys
function staticAtoms(orig_atoms::ExtXYZ.Atoms)
    new_pos  = [SVector{3}(pos) for pos in position(orig_atoms)]
    new_vels = [SVector{3}(pos) for pos in velocity(orig_atoms)]

    new_atoms_data = Dict(:atomic_symbol => atomic_symbol(orig_atoms),
                          :atomic_number => atomic_number(orig_atoms),
                          :atomic_mass   => atomic_mass(orig_atoms), 
                          :position     => new_pos,
                          :velocity      => new_vels)
    
    ExtXYZ.Atoms(NamedTuple(new_atoms_data),orig_atoms.system_data)
end
