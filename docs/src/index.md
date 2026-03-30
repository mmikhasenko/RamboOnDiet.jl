# RemboOnDiet

`RemboOnDiet.jl` implements a RAMBO-on-diet style phase-space generator for relativistic
multi-particle final states. It uses [`FourVectors.jl`](https://github.com/mmikhasenko/FourVectors.jl)
and `LorentzVectorBase.jl` for four-vector storage and operations, and it exposes a single
`PhaseSpaceGenerator` object that can be sampled with `rand(generator)`.

The package is organized around three ideas:

- a fixed-size map from `3n - 4` unit-hypercube variables to an `n`-body event,
- exact four-momentum conservation by construction through sequential two-body decays,
- a phase-space weight that is constant in the massless case and becomes the corrective
  Jacobian for massive final states.

```@contents
Pages = ["algorithm.md", "implementation.md", "generated/three-body.md"]
Depth = 2
```

## Guide

- [Algorithm overview](algorithm.md)
- [Implementation details](implementation.md)
- [Three-body generation and Dalitz processing](generated/three-body.md)
