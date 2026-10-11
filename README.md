# RamboOnDiet.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/RamboOnDiet.jl/)
[![CI](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml)
[![Docs workflow](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml)

`RamboOnDiet.jl` is a Julia package for generating relativistic `N`-body
phase-space points with a RAMBO-on-diet style sampler, with emphasis on
numerically validated phase-space generation and Dalitz-plot studies.

The original version of the package was implemented using Codex/GPT-5.4,
based on the [RAMBO on diet](https://arxiv.org/abs/1308.2922) research paper,
validated using `ThreeBodyDecays.jl`, and reviewed by a human.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/mmikhasenko/RamboOnDiet.jl")
```

Dependencies such as [`FourVectors.jl`](https://github.com/JuliaHEP/FourVectors.jl)
and [`ThreeBodyDecays.jl`](https://github.com/mmikhasenko/ThreeBodyDecays.jl)
are registered in General and resolved automatically.

## Docs

- [Documentation website](https://mmikhasenko.github.io/RamboOnDiet.jl/)
- [Documentation overview](docs/src/index.md)
- [Three-body tutorial](docs/src/generated/three-body.md)

To build the rendered docs locally, run:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

## Quick Start

```julia
using RamboOnDiet
using FourVectors

masses = [0.93827208816, 0.493677, 0.13957039]
generator = PhaseSpaceGenerator(masses, 2.28646)
point = rand(generator)

weight = phase_space_weight(point)
total = total_momentum(point.momenta)
```

## Aligned phase space

For `n >= 3`, fix the overall orientation with particle 1 along -z and particle 2
in the xz plane with nonnegative x, in the total rest frame:

```julia
point = generate_aligned_point(generator)
point = generate_aligned_from_unit_hypercube([0.4, 0.6], generator) # three bodies
```

The reduced map uses `3n - 7` unit-cube coordinates: `n - 2` mass coordinates,
then the polar coordinate for particle 2, then polar/azimuthal pairs for subsequent
splittings. Weights include the integrated overall orientation. See the
[coordinate maps and aligned-vector conventions](docs/src/algorithm.md#aligned-four-vectors).
