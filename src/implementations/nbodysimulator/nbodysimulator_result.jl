# -----------------------------------------------------------------------------
# Implementation of Atomistic MolecularDynamicsResult API
# -----------------------------------------------------------------------------

"""
    NBSResult <: MolecularDynamicsResult

The result generated from running an `NBSimulator`.

**Field descriptions**
- `result::SimulationResult` the standard simulation result from `NBodySimulator`
"""
struct NBSResult <: MolecularDynamicsResult
    result::SimulationResult
end

get_time_range(result::NBSResult) = result.result.solution.t * TIME_UNIT
get_bounding_box(result::NBSResult) = get_bounding_box(result.result.simulation.boundary_conditions)
get_boundary_conditions(result::NBSResult) = get_boundary_conditions(result.result.simulation.boundary_conditions)

reference_temperature(result::NBSResult) = reference_temperature(result.result.simulation.thermostat)
reference_temperature(thermostat::Thermostat) = thermostat.T * TEMPERATURE_UNIT
reference_temperature(::NullThermostat) = missing

function get_positions(result::NBSResult, t::Integer)
    positions = get_position(result.result, result.result.solution.t[t])
    [bound_position(SVector{3}(p), result.result.simulation.boundary_conditions) for p ∈ eachcol(positions)] * LENGTH_UNIT
end
function get_velocities(result::NBSResult, t::Integer)
    velocities = get_velocity(result.result, result.result.solution.t[t])
    [SVector{3}(v) for v ∈ eachcol(velocities)] * VELOCITY_UNIT
end
function get_particles(result::NBSResult, t::Integer)
    [Atom(b.symbol, p, v; b.data...) for (b, p, v) ∈ zip(result.result.simulation.system.bodies, get_positions(result, t), get_velocities(result, t))]
end

temperature(result::NBSResult, t::Integer) = NBodySimulator.temperature(result.result, result.result.solution.t[t]) * TEMPERATURE_UNIT
kinetic_energy(result::NBSResult, t::Integer) = NBodySimulator.kinetic_energy(result.result, result.result.solution.t[t]) * ENERGY_UNIT
function potential_energy(result::NBSResult, t::Integer)
    potentials = result.result.simulation.system.potentials
    # https://github.com/SciML/NBodySimulator.jl/issues/44
    if :custom ∈ keys(potentials)
        InteratomicPotentials.potential_energy(get_system(result, t), potentials[:custom].potential) * ENERGY_UNIT
    else
        NBodySimulator.potential_energy(result.result, result.result.solution.t[t]) * ENERGY_UNIT
    end
end

function rdf(result::NBSResult, sample_fraction::Float64 = 1.0)
    @assert 0 < sample_fraction ≤ 1
    sr = result.result
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    trange = sr.solution.t[end-floor(Int, length(sr.solution.t) * sample_fraction)+1:end]

    maxbin = 1000
    radius = 0.5 * pbc.L
    dr = radius / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for i ∈ 1:n
            ri = SVector{3}(cc[:, i])
            for j ∈ i+1:n
                rj = SVector{3}(cc[:, j])
                (rij, r, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)
                if r < radius
                    bin = ceil(Int, r / dr)
                    hist[bin] += 2
                end
            end
        end
    end

    c = 4 / 3 * π * n / pbc.L^3
    gr = zeros(maxbin)
    rs = zeros(maxbin)
    tlen = length(trange)
    for bin ∈ 1:maxbin
        rlower = (bin - 1) * dr
        rupper = rlower + dr
        nideal = c * (rupper^3 - rlower^3)
        gr[bin] = (hist[bin] / (tlen * n)) / nideal
        rs[bin] = rlower + dr / 2
    end

    rs, gr
end
