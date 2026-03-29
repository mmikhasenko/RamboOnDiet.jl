module PhaseSpaceRembo

using LinearAlgebra
using Random
using StaticArrays
using FourVectors
using LorentzVectorBase

export PhaseSpacePoint
export PhaseSpaceGenerator
export generate_point
export generate_momenta
export generate_from_unit_hypercube
export required_random_numbers
export phase_space_weight
export total_momentum
export mandelstam
export invariant_masses
export breakup_momentum
export kallen
export massless_phase_space_volume

include("kinematics.jl")
include("solver.jl")
include("kinematics_object.jl")
include("generator.jl")

end
