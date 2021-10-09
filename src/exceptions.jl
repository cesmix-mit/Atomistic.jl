# Exceptions

struct UnimplementedError <: Exception
    func::Symbol
    type::Type
    UnimplementedError(func::Symbol, any) = new(func, typeof(any))
end

Base.showerror(io::IO, e::UnimplementedError) = print(io, e.type, " doesn't implement ", e.func)