# Unit tests for unit_convention.jl

import Unitful.𝐌, Unitful.𝐋, Unitful.𝐓, Unitful.𝚯

@testset "unit_convention.jl" begin
    @test Atomistic.MASS_UNIT == aunit(𝐌)
    @test Atomistic.LENGTH_UNIT == aunit(𝐋)
    @test Atomistic.ENERGY_UNIT == aunit(𝐌 * 𝐋^2 / 𝐓^2)
    @test Atomistic.TIME_UNIT == aunit(𝐓)
    @test Atomistic.VELOCITY_UNIT == aunit(𝐋 / 𝐓)
    @test Atomistic.TEMPERATURE_UNIT == aunit(𝚯)
    @test Atomistic.FORCE_UNIT == aunit(𝐌 * 𝐋 / 𝐓^2)

    @test Atomistic.MASS_TYPE == typeof(1.0aunit(𝐌))
    @test Atomistic.LENGTH_TYPE == typeof(1.0aunit(𝐋))
    @test Atomistic.ENERGY_TYPE == typeof(1.0aunit(𝐌 * 𝐋^2 / 𝐓^2))
    @test Atomistic.TIME_TYPE == typeof(1.0aunit(𝐓))
    @test Atomistic.VELOCITY_TYPE == typeof(1.0aunit(𝐋 / 𝐓))
    @test Atomistic.TEMPERATURE_TYPE == typeof(1.0aunit(𝚯))
    @test Atomistic.FORCE_TYPE == typeof(1.0aunit(𝐌 * 𝐋 / 𝐓^2))
end
