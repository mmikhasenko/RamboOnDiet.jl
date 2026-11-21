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

export PhaseSpacePoint
export PhaseSpaceGenerator
export generate_point
export generate_from_unit_hypercube
export generate_aligned_point
export generate_aligned_from_unit_hypercube
export phase_space_weight

include("kinematics.jl")
include("solver.jl")
include("kinematics_object.jl")
include("generator.jl")

end
