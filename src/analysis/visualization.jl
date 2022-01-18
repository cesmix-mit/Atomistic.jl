# Functions for exporting molecular dynamics simulation results to visualization tools

"""
    write_nbs_animation(result::NBSResult, filename::String)

Animate an `NBSResult` and store the result in a .gif file.
"""
function write_nbs_animation(result::NBSResult, filename::String)
    animate(result.result, filename)
end

"""
    write_ase_trajectory(result::MolecularDynamicsResult, element::DFTK.Element, lattice, filename::String)

Write the trajectory of a `MolecularDynamicsResult` to a .traj file.

The file can be visualized by running `ase gui <filename>` on the command line.
"""
function write_ase_trajectory(result::MolecularDynamicsResult, filename::String)
    ase = pyimport_e("ase")
    if ispynull(ase)
        @warn "ASE is not installed, skipping write_trajectory"
    else
        filename = abspath(expanduser(filename))
        if isfile(filename)
            @info "overwriting trajectory at" filename
            rm(filename)
        else
            @info "writing trajectory to" filename
        end
        traj = pyimport("ase.io.trajectory").Trajectory(filename, "a")
        for t ∈ 1:length(get_time_range(result))
            lattice, atoms = parse_system(get_system(result, t))
            traj.write(ase_atoms(lattice, atoms))
        end
        traj.close()
    end
end
