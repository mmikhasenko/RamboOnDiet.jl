"""
    RamboOnDiet

Relativistic `n`-body phase-space sampler using the RAMBO-on-diet map.
Daughter masses and the total energy are in one common energy unit.

Load the package and build a [`PhaseSpaceGenerator`](@ref). The masses are
`generator.masses`; draw an event from that object next.

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
generator.masses
```

## Related

- [`PhaseSpaceGenerator`](@ref) — set the masses and the total four-momentum
"""
module RamboOnDiet

using LinearAlgebra
using Random
using StaticArrays
using FourVectors
using LorentzVectorBase

# Kept: a user builds a generator, draws or maps an event, then reads its weight.
export PhaseSpacePoint # kept: event the samplers return; read point.momenta and point.weight
export PhaseSpaceGenerator # kept: entry point; masses and the total four-momentum
export generate_point # kept: draw one event; rand(generator) calls this
export generate_from_unit_hypercube # kept: map 3n-4 unit-cube coordinates to one event
export generate_aligned_point # kept: draw one event with the overall rotation fixed
export generate_aligned_from_unit_hypercube # kept: map the reduced aligned coordinates
export phase_space_weight # kept: read the Monte Carlo weight stored on an event

# Unexported names stay defined for the tests and for the sampler.
#   generate_momenta — generate_point on PhaseSpaceGenerator(masses, total)
#   required_random_numbers — generator.random_numbers
#   required_aligned_random_numbers — 0 for two bodies, otherwise 3n-7
#   total_momentum — sum of point.momenta
#   mandelstam, invariant_masses — mass2 of a sum of daughters
#   resample_mass_zeros!, threshold_configuration, solve_mass_parameter,
#   breakup_momentum, two_body_density, massless_phase_space_volume,
#   boost, unit_direction, kallen, spatial, fourvector, add_fourvectors, Vec3
#     — steps inside the map

include("kinematics.jl")
include("solver.jl")
include("kinematics_object.jl")
include("generator.jl")

end
