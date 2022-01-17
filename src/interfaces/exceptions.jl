# Exceptions

"""
    UnimplementedError
Exception thrown in default implementation of API to indicate that an implementator did not provide an implementation of a particular API function.
"""
struct UnimplementedError <: Exception
    func::Symbol
    type::Type
    UnimplementedError(func::Symbol, any) = new(func, typeof(any))
end

Base.showerror(io::IO, e::UnimplementedError) = print(io, e.type, " doesn't implement ", e.func)