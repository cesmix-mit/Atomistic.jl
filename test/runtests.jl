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
        include("analysis/visualization.jl")
    end
end
