# RamboOnDiet

`RamboOnDiet.jl` implements a RAMBO-on-diet style phase-space generator for relativistic
multi-particle final states. It uses [`FourVectors.jl`](https://github.com/JuliaHEP/FourVectors.jl)
and [`LorentzVectorBase.jl`](https://github.com/JuliaHEP/LorentzVectorBase.jl) for four-vector storage and operations, and it exposes a single
`PhaseSpaceGenerator` object that can be sampled with `rand(generator)`.

The package is organized around three ideas:

- a fixed-size map from `3n - 4` unit-hypercube variables to an `n`-body event,
- exact four-momentum conservation by construction through sequential two-body decays,
- a phase-space weight that is constant in the massless case and becomes the corrective
  Jacobian for massive final states.

The implementation always uses the same main sampler. The 2-body and 3-body cases are
validation anchors because their expected geometry is especially transparent, not because
they use a separate generator.

## Quick Start

```julia
using RamboOnDiet

generator = PhaseSpaceGenerator([0.1, 0.2, 0.3, 0.4], 5.0)
point = rand(generator)

point.momenta
phase_space_weight(point)
```

## Main API

The package exports eight names:

| Name | Purpose |
|:--|:--|
| `PhaseSpaceGenerator` | Configure masses and total energy or four-momentum |
| `PhaseSpacePoint` | Store the generated four-vectors and event weight |
| `generate_point` | Sample an ordinary point; also available as `rand(generator)` |
| `decay_two_body` | Construct a two-body point at explicit rest-frame angles |
| `generate_from_unit_hypercube` | Map `3n - 4` unit-cube inputs to a point |
| `generate_aligned_point` | Sample a point with fixed overall orientation |
| `generate_aligned_from_unit_hypercube` | Map zero (two-body) or `3n - 7` inputs to an aligned point |
| `phase_space_weight` | Read the phase-space integration weight |

Use `PhaseSpaceGenerator(masses, sqrt_s)` for a parent at rest, or
`PhaseSpaceGenerator(masses, total::FourVector)` for a moving parent.
Read momenta through `point[i]` or `point.momenta`. The ordinary coordinate count
is `generator.random_numbers`; the aligned count is zero for two bodies and `3 * length(generator) - 7` otherwise.

Kinematic analysis uses the existing four-vector operations:

```julia
using FourVectors

sum(point.momenta)          # total four-momentum
mass2(point[2] + point[3])   # pair invariant s23 = (p2 + p3)^2
```

```@docs
decay_two_body
```

### Migrating from the larger export list

Coordinate-count helpers, kinematic helpers, and `generate_momenta` are no longer
exported. Use the expressions above for counts, totals, and pair invariants; replace
`generate_momenta(masses, total)` with `rand(PhaseSpaceGenerator(masses, total))`
(or `rand(rng, generator)` with an explicit RNG). The old helper implementations
remain accessible with the `RamboOnDiet.` prefix for compatibility, but are internal
and are not part of the supported public API. In particular, the old `mandelstam`
helper computed a pair invariant mass squared, now written directly as
`mass2(point[i] + point[j])`.

```@contents
Pages = ["algorithm.md", "implementation.md", "generated/three-body.md", "generated/b-decay.md"]
Depth = 2
```

## Guide

- [Algorithm overview](algorithm.md)
- [Implementation details](implementation.md)
- [Three-body generation and Dalitz processing](generated/three-body.md)
- [B decay cascade and Arrow export](generated/b-decay.md)
