# [WIP] Atomistic.jl

Package that provides an integrated workflow for molecular dyanmics simulations. Defines an API for molecular dynamics (MD) simulations that is compatible with the interatomic potential interface defined by [Potentials.jl](https://github.com/cesmix-mit/Potentials.jl) and the atomic configuration interface defined by [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl).

Developed as part of the [CESMIX](https://computing.mit.edu/cesmix) Julia package suite. See also [Potentials.jl](https://github.com/cesmix-mit/Potentials.jl), [PotentialLearning.jl](https://github.com/cesmix-mit/PotentialLearning.jl), and [PotentialUQ.jl](https://github.com/cesmix-mit/PotentialUQ.jl).

If you want to integrate an existing MD code with the Atomistic API, see [Implementing the Atomistic API](@ref). If you want to use a code that is already integrated with Atomistic to run MD simulations see [Using Atomistic-Compatible Packages](@ref).
