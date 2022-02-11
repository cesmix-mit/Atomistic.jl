# Unit tests for unit_convention.jl

@testset "unit_convention.jl" begin
    @test Atomistic.MASS_UNIT == aunit(u"kg")
    @test Atomistic.LENGTH_UNIT == aunit(u"m")
    @test Atomistic.ENERGY_UNIT == aunit(u"J")
    @test Atomistic.TIME_UNIT == aunit(u"s")
    @test Atomistic.VELOCITY_UNIT == aunit(u"m/s")
    @test Atomistic.TEMPERATURE_UNIT == aunit(u"K")
end
