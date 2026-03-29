# RemboOnDiet.jl

`RemboOnDiet.jl` is a Julia package for generating relativistic `N`-body
phase-space points with a RAMBO-on-diet style sampler.

## Docs

- [Documentation overview](docs/src/index.md)
- [Three-body tutorial](docs/src/generated/three-body.md)

To build the rendered docs locally, run:

```bash
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

## Quick Start

```julia
using RemboOnDiet

masses = [0.93827208816, 0.493677, 0.13957039]
generator = PhaseSpaceGenerator(masses, 2.28646)
point = rand(generator)

weight = phase_space_weight(point)
total = total_momentum(point.momenta)
```
