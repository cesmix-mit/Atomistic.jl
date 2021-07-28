function write_trajectory(result::MolecularDynamicsResult, element::DFTK.Element, lattice::AbstractArray{Quantity, 2}, filename::String)
    if isfile(filename)
        println("replacing $(filename)")
        rm(filename)
    end
    traj = pyimport("ase.io.trajectory").Trajectory(filename, "a")
    for t ∈ 1:length(get_time_range(result))
        state = get_bodies(result, t)
        traj.write(ASEAtoms(state, element, lattice).atoms)
    end
    traj.close()
end
