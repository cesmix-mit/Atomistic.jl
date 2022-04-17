# Unit tests for api/exceptions.jl

@testset "api/exceptions.jl" begin
    target = UnimplementedError(:func, 42)

    @test target isa Exception
    @test target isa UnimplementedError{Int64}

    message = sprint(showerror, target)
    @test contains(message, "func")
    @test contains(message, "Int64")
end
