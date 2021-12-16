# Tests for Atomistic.jl

using Atomistic
using PeriodicTable
using Test
using Unitful
using UnitfulAtomic

@testset "Atomistic.jl" begin
    @testset "Unit Tests" begin
        include("unit_tests_ab_integration.jl")
        include("unit_tests_nbs_integration.jl")
        include("unit_tests_dftk_integration.jl")
    end
    @testset "Integration Tests" begin
        include("integration_test_ab_nbs.jl")
        include("integration_test_ab_nbs_ip.jl")
        include("integration_test_ab_nbs_dftk.jl")
    end
end
