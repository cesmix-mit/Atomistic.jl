# Implementation of Atomistic API with Molly.jl

#=
The contents of molly_base.jl have largely been superceded by the native
Molly.System(::AbstractSystem). However the new molly_temporary_utils.jl file is needed to
ensure that the ExtXYZ.Atoms objects (returned by AtomsIO) has static arrays for pos and vels

So one of two longer-term solutions are needed
1. Modifying the validate_coords function in Molly to be able to handle Vector{Vector}
   instead of just Vector{SArray}
2. ExtXYZ.Atoms should use StaticArrays for positions and velocities
=#
#include("molly/molly_base.jl")
include("molly/molly_temporary_utils.jl")
#include("molly/molly_simulator.jl")
#include("molly/molly_result.jl")
include("molly/interpot_wrapper.jl")
