# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------
const E_DIM = dimension(u"eV")
const L_DIM = dimension(u"Å")

const DEF_EUNIT = u"eV"
const DEF_LUNIT = u"Å"

const NOUNIT_T = typeof(NoUnits)

struct InteratomicPotentialInter{P<:AbstractPotential}
    potential::P
    energy_units::Unitful.Unitlike 
    length_units::Unitful.Unitlike

    # internal constructor, ensuring energy units and length units have correct dimensions
    ( InteratomicPotentialInter(pot::AbstractPotential, 
                                eu::Union{NOUNIT_T, Unitful.Units{UE,E_DIM,nothing}} = DEF_EUNIT,
                                lu::Union{NOUNIT_T, Unitful.Units{UL,L_DIM,nothing}} = DEF_LUNIT)
                                where {UE,UL} = new{typeof(pot)}(pot,eu,lu) )
end

function Molly.forces(inter::InteratomicPotentialInter, 
                      sys::AbstractSystem,
                      neighbors = nothing;
                      n_threads = Threads.nthreads())
    forces = InteratomicPotentials.force(sys,inter.potential)

    # initial profiling didn't show huge performance hit from unit conversion 
    # but! that may be because other parts of IP.jl are very slow, e.g. neighbor list construction
    if eltype(forces[1]) <: Unitful.Quantity 
        if inter.energy_units != NoUnits
            forces = [uconvert.(inter.energy_units/inter.length_units, fi)
                     for fi in forces]
        else
            forces = [ustrip.(fi)
                     for fi in forces]
        end
    elseif eltype(forces[1]) <: Real && inter.energy_units != NoUnits
        forces = forces * inter.energy_units/inter.length_units
    end

    forces
end

function Molly.potential_energy(inter::InteratomicPotentialInter, 
                                sys::AbstractSystem,
                                neighbors = nothing;
                                n_threads = Threads.nthreads())

    energy = InteratomicPotentials.potential_energy(sys,inter.potential)

    if typeof(energy) <: Unitful.Quantity 
        if inter.energy_units != NoUnits
            energy = uconvert(inter.energy_units, energy)
        else
            energy = ustrip(energy)
        end
    elseif typeof(energy) <: Real && inter.energy_units != NoUnits
        energy = energy * inter.energy_units
    end

    energy
end

