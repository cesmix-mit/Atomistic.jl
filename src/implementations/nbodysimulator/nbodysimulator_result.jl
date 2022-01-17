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

function get_system(result::NBSResult, t::Integer = 0)
    sr = result.result
    simulation = sr.simulation
    system = simulation.system
    time = sr.solution.t[t > 0 ? t : end]
    positions = get_position(sr, time)
    velocities = get_velocity(sr, time)
    masses = get_masses(system)
    symbols = get_atomic_symbols(system)
    numbers = get_atomic_numbers(system)
    particles = [ElementMassBody(SVector{3}(positions[:, i]), SVector{3}(velocities[:, i]), masses[i], symbols[i], numbers[i]) for i ∈ 1:length(system.bodies)]
    DynamicSystem(FlexibleSystem(particles, simulation.boundary_conditions), time * TIME_UNIT)
end

function get_time_range(result::NBSResult)
    result.result.solution.t * TIME_UNIT
end

function temperature(result::NBSResult, t::Integer = 0)
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.temperature(result.result, time) * TEMPERATURE_UNIT
end

function reference_temperature(result::NBSResult)
    thermostat = result.result.simulation.thermostat
    thermostat isa NullThermostat ? missing : thermostat.T * TEMPERATURE_UNIT
end

function kinetic_energy(result::NBSResult, t::Integer = 0)
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.kinetic_energy(result.result, time) * ENERGY_UNIT
end

function potential_energy(result::NBSResult, t::Integer = 0)
    potentials = result.result.simulation.system.potentials
    # https://github.com/SciML/NBodySimulator.jl/issues/44
    if :custom ∈ keys(potentials)
        return InteratomicPotentials.potential_energy(get_system(result, t), potentials[:custom].potential) * ENERGY_UNIT
    end
    time = result.result.solution.t[t > 0 ? t : end]
    NBodySimulator.potential_energy(result.result, time) * ENERGY_UNIT
end

function rdf(result::NBSResult, sample_fraction::Float64 = 1.0)
    @assert 0 < sample_fraction ≤ 1
    sr = result.result
    n = length(sr.simulation.system.bodies)
    pbc = sr.simulation.boundary_conditions
    trange = sr.solution.t[end-floor(Int, length(sr.solution.t) * sample_fraction)+1:end]

    maxbin = 1000
    dr = pbc.L / 2 / maxbin
    hist = zeros(maxbin)
    for t ∈ trange
        cc = get_position(sr, t)
        for i ∈ 1:n
            ri = SVector{3}(cc[:, i])
            for j ∈ i+1:n
                rj = SVector{3}(cc[:, j])
                # TODO: potential for major performance improvements here
                (rij, r, r2) = NBodySimulator.get_interparticle_distance(ri, rj, pbc)
                if r2 < (0.5 * pbc.L)^2
                    bin = ceil(Int, r / dr)
                    if bin > 1 && bin <= maxbin
                        hist[bin] += 2
                    end
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
