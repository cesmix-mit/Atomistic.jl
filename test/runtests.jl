# Tests for Atomistic.jl

using Atomistic
using AtomsBase
using InteratomicPotentials
using PeriodicTable
using Test
using Unitful
using UnitfulAtomic

@testset "Atomistic.jl" begin
    @testset "Unit Tests" begin
        include("unit/unit_convention.jl")
        @testset "Atomistic API Unit Tests" begin
            include("unit/api/exceptions.jl")
            include("unit/api/md_result.jl")
            include("unit/api/md_simulator.jl")
        end
        @testset "Integrations Unit Tests" begin
            include("unit/integrations/atomsbase_integration.jl")
            include("unit/integrations/dftk_integration.jl")
        end
        @testset "Implementations Unit Tests" begin
            @testset "NBodySimulator Implementation" begin
                include("unit/implementations/nbodysimulator/nbodysimulator_base.jl")
                include("unit/implementations/nbodysimulator/nbodysimulator_simulator.jl")
                include("unit/implementations/nbodysimulator/nbodysimulator_result.jl")
            end
        end
        @testset "Analysis Unit Tests" begin
            include("unit/analysis/plotting.jl")
            include("unit/analysis/visualization.jl")
        end
    end
    @testset "Integration Tests" begin
        include("integration/atomsbase_nbodysimulator.jl")
        include("integration/atomsbase_nbodysimulator_interatomicpotentials.jl")
        include("integration/atomsbase_nbodysimulator_dftk.jl")
    end
end
