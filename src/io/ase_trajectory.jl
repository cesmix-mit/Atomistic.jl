function write_trajectory(result::MolecularDynamicsResult, box_size::Quantity, element::DFTK.Element, lattice::AbstractArray{Quantity, 2}, filename::String)
    if isfile(filename)
        println("replacing $(filename)")
        rm(filename)
    end
    traj = pyimport("ase.io.trajectory").Trajectory(filename, "a")
    for t ∈ 1:length(get_time_range(result))
        bodies = get_bodies(result, t)
        atoms = ase_atoms(element, bodies, box_size, lattice)
        traj.write(atoms)
    end
    traj.close()
end
