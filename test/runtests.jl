# Tests for Atomistic.jl

using Atomistic

using AtomsBase
using DFTK
using InteratomicPotentials
using Molly
using NBodySimulator
using StaticArrays
using Test
using Unitful

Base.:(==)(a1::AtomsBase.Atom, a2::AtomsBase.Atom) = propertynames(a1) == propertynames(a2) && all(getproperty(a1, p) == getproperty(a2, p) for p ∈ propertynames(a1))
Base.:(==)(b1::Atomistic.ElementMassBody, b2::Atomistic.ElementMassBody) = all(getproperty(b1, p) == getproperty(b2, p) for p ∈ propertynames(b1))
Base.:(==)(a1::Atomistic.AugmentedAtomData, a2::Atomistic.AugmentedAtomData) = propertynames(a1) == propertynames(a2) && all(getproperty(a1, p) == getproperty(a2, p) for p ∈ propertynames(a1))

@testset "Atomistic.jl" begin
    include("unit_convention.jl")
    @testset "Atomistic API Unit Tests" begin
        include("api/exceptions.jl")
        include("api/md_result.jl")
        include("api/md_simulator.jl")
    end
    @testset "Integrations Unit Tests" begin
        include("integrations/initialization.jl")
        include("integrations/dftk_integration.jl")
    end
    @testset "Implementations Unit Tests" begin
        @testset "NBodySimulator Implementation" begin
            include("implementations/nbodysimulator/nbodysimulator_base.jl")
            include("implementations/nbodysimulator/nbodysimulator_simulator.jl")
            include("implementations/nbodysimulator/nbodysimulator_result.jl")
        end
        @testset "Molly Implementation" begin
            include("implementations/molly/molly_base.jl")
            include("implementations/molly/molly_simulator.jl")
            include("implementations/molly/molly_result.jl")
        end
    end
    @testset "Analysis Unit Tests" begin
        include("analysis/plotting.jl")
    end
end
