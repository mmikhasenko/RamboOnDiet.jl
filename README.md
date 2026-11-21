# RamboOnDiet.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/RamboOnDiet.jl/)
[![CI](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml)
[![Docs workflow](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml)

`RamboOnDiet.jl` samples relativistic `n`-body phase space with a RAMBO-on-diet map.
Daughter masses and the total energy use one common energy unit. The map follows
[RAMBO on diet](https://arxiv.org/abs/1308.2922).

Build a `PhaseSpaceGenerator`, draw an event with `rand` or `generate_point`, and read
`phase_space_weight` before using the momenta. Chosen coordinates go through
`generate_from_unit_hypercube`. Observables invariant under the overall rotation
use `generate_aligned_point` or `generate_aligned_from_unit_hypercube`.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/mmikhasenko/RamboOnDiet.jl")
```

[`FourVectors.jl`](https://github.com/JuliaHEP/FourVectors.jl) is a dependency.
[`ThreeBodyDecays.jl`](https://github.com/mmikhasenko/ThreeBodyDecays.jl) is used
in the tests and tutorials. Both are registered in General.

## Sampling

```julia
using RamboOnDiet
using FourVectors

masses = [0.93827208816, 0.493677, 0.13957039]
generator = PhaseSpaceGenerator(masses, 2.28646)
point = rand(generator)

weight = phase_space_weight(point)
total = sum(point.momenta)
```

`generator.masses` lists the daughters. `generator.invariant_mass` and
`generator.threshold` are the parent mass and the sum of the daughter masses.
`generator.random_numbers` is `3n - 4`, the length of the ordinary unit-cube vector.
`point.momenta` and `point[i]` are the daughter four-momenta.

Pass a `FourVector` when the parent is moving:

```julia
using RamboOnDiet, FourVectors
masses = [0.93827208816, 0.493677, 0.13957039]
total = FourVector(0.4, -0.3, 1.2; E = sqrt(2.28646^2 + 0.4^2 + 0.3^2 + 1.2^2))
point = generate_point(PhaseSpaceGenerator(masses, total))
point.momenta
```

A fixed coordinate vector uses the same generator. Above threshold, the first
`n - 2` entries lie in `(0, 1)` and the remaining angles lie in `[0, 1]`.

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
point = generate_from_unit_hypercube([0.3, 0.2, 0.7, 0.4, 0.9], generator)
point.weight
```

A pair invariant is `mass2(point[i] + point[j])` after `using FourVectors`.

## Aligned phase space

Two-body alignment needs no coordinates: particle 1 points along +z and particle 2
along -z in the center-of-mass frame.

```julia
using RamboOnDiet
point = generate_aligned_from_unit_hypercube(Float64[], PhaseSpaceGenerator([0.3, 0.5], 2.0))
```

For `n >= 3`, particle 1 points along -z and particle 2 lies in the xz plane with
nonnegative x, in the total rest frame. The reduced map uses `3n - 7` unit-cube
coordinates: `n - 2` mass coordinates, then the polar coordinate for particle 2,
then polar/azimuthal pairs for later splittings. The weight includes the integrated
overall orientation.

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.93827208816, 0.493677, 0.13957039], 2.28646)
point = generate_aligned_point(generator)
point = generate_aligned_from_unit_hypercube([0.4, 0.6], generator) # three bodies
```

The coordinate tables are in the
[algorithm notes](docs/src/algorithm.md#aligned-four-vectors).
Random samplers redraw an exact-zero mass input. Deterministic maps throw
`DomainError` for a mass endpoint. Angular coordinates may include both endpoints.

## Documentation

- [Documentation website](https://mmikhasenko.github.io/RamboOnDiet.jl/)
- [Manual overview](docs/src/index.md)
- [Three-body tutorial](docs/src/generated/three-body.md)
- [B decay cascade and Arrow export](docs/src/generated/b-decay.md)

Build the rendered manual locally with:

```bash
julia --project=docs docs/generate.jl
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`. After editing a tutorial `.jl` file, regenerate
and commit its Markdown in `docs/src/generated/`. CI runs only `docs/make.jl`.
