# Integrations with InteratomicPotentials.jl

function InteratomicPotentials.force(r::MassBodies, p::ArbitraryPotential)
    force([Atom(b.m, Vector(mod.(b.r, austrip(r.box_size))), Vector(b.v)) for b in r.bodies], p)
end