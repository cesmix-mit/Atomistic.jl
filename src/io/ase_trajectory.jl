function write_trajectory(result::NBodySimulator.SimulationResult, box_size::Quantity, element::DFTK.Element, lattice::AbstractArray{Quantity, 2}, filename::String)
    if isfile(filename)
        println("replacing $(filename)")
        rm(filename)
    end
    traj = pyimport("ase.io.trajectory").Trajectory(filename, "a")
    for t ∈ 1:length(result.solution.t)
        bodies = extract_bodies(result, t)
        atoms = dftk_atoms(element, bodies, box_size)
        asa = ase_atoms(austrip.(lattice), atoms)
        traj.write(asa)
    end
    traj.close()
end
