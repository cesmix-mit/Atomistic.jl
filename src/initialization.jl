# Convenience functions for initializing AtomsBase systems.

# See NBodySimulator.generate_bodies_in_cell_nodes
function generate_atoms_in_cubic_cell(n::Integer, symbol::Union{Integer,AbstractString,Symbol}, L::Unitful.Length, reference_temp::Unitful.Temperature)
    average_velocity = √(u"k" * reference_temp / elements[symbol].atomic_mass)
    velocities = average_velocity * randn(Float64, (3, n))
    particles = Vector{AtomsBase.Atom}(undef, n)

    count = 1
    dL = L / (ceil(n^(1 / 3)))
    for x = dL/2:dL:L, y = dL/2:dL:L, z = dL/2:dL:L
        if count > n
            break
        end
        particles[count] = AtomsBase.Atom(symbol, SVector{3}(x, y, z), SVector{3}(velocities[:, count]))
        count += 1
    end

    bounding_box = SVector{3}(SVector{3}.(eachrow(L * I(3))))
    periodic_system(particles, bounding_box)
end
