# [WIP] Atomistic.jl

Package that provides an integrated workflow for molecular dyanmics simulations. Defines an API for molecular dynamics (MD) simulations that is compatible with the interatomic potential interface defined by [Interatomicotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl) and the atomic configuration interface defined by [AtomsBase.jl](https://github.com/JuliaMolSim/AtomsBase.jl).

Developed as part of the [CESMIX](https://computing.mit.edu/cesmix/) Julia package suite. See also [InteratomicPotentials.jl](https://github.com/cesmix-mit/InteratomicPotentials.jl), [PotentialLearning.jl](https://github.com/cesmix-mit/PotentialLearning.jl), and [PotentialUQ.jl](https://github.com/cesmix-mit/PotentialUQ.jl).

## Conventions

The unit convention throughout the package and other packages in the CESMIX Julia package suite is to assume all unspecified units to be atomic units as defined in the [UnitfulAtomic](https://github.com/sostock/UnitfulAtomic.jl) package. All exposed interfaces should allow for numeric or unitful input. For clarity's sake, it is _strongly recommended_ that user code utilize unitful wherever possible. Internally, Atomistic will automatically convert these quantities to be compatible without any significant performance penalty.

## Next Steps

If you want to integrate an existing MD code with the Atomistic API, see [Implementing the Atomistic API](@ref). If you want to use a code that is already integrated with Atomistic to run MD simulations, see [Using Atomistic-Compatible Packages](@ref). If you want to see the full API reference (for Atomistic API and other types and functions exported by the package), see [API Reference](@ref).
