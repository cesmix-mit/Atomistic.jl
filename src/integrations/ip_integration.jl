# Integrations with InteratomicPotentials.jl

function InteratomicPotentials.force(r::MassBodies, p::ArbitraryPotential)
    force([Atom(b.m, Vector(b.r), Vector(b.v)) for b in r.bodies], p)
end