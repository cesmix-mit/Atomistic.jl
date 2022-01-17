# Exceptions

"""
    UnimplementedError
Exception thrown in default implementation of API to indicate that an implementator did not provide an implementation of a particular API function.
"""
struct UnimplementedError{T} <: Exception
    func::Symbol
end
UnimplementedError(func::Symbol, obj) = UnimplementedError{typeof(obj)}(func)

Base.showerror(io::IO, e::UnimplementedError{T}) where {T} = print(io, T, " doesn't implement ", e.func)