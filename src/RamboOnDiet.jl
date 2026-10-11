module RamboOnDiet

using LinearAlgebra
using Random
using StaticArrays
using FourVectors
using LorentzVectorBase

export PhaseSpacePoint
export PhaseSpaceGenerator
export generate_point
export decay_two_body
export generate_from_unit_hypercube
export generate_aligned_point
export generate_aligned_from_unit_hypercube
export phase_space_weight

include("kinematics.jl")
include("solver.jl")
include("kinematics_object.jl")
include("generator.jl")

end
