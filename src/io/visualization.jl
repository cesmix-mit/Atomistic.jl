function write_nbs_animation(result::NBSResult, filename::String)
    animate(result.simulation_result, filename)
end

function write_ase_trajectory(result::MolecularDynamicsResult, element::DFTK.Element, lattice::AbstractArray{Quantity, 2}, filename::String)
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
            state = get_bodies(result, t)
            traj.write(ASEAtoms(state, element, lattice).atoms)
        end
        traj.close()
    end
end
