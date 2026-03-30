# RemboOnDiet.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/RamboOnDiet.jl/)
[![CI](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/ci.yml)
[![Docs workflow](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml/badge.svg)](https://github.com/mmikhasenko/RamboOnDiet.jl/actions/workflows/docs.yml)

`RemboOnDiet.jl` is a Julia package for generating relativistic `N`-body
phase-space points with a RAMBO-on-diet style sampler, with emphasis on
numerically validated phase-space generation and Dalitz-plot studies.

The original version of the package was implemented using Codex/GPT-5.4,
based on the [RAMBO on diet](https://arxiv.org/abs/1308.2922) research paper,
validated using `ThreeBodyDecays.jl`, and reviewed by a human.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/mmikhasenko/FourVectors.jl")
Pkg.add(url = "https://github.com/mmikhasenko/RamboOnDiet.jl")
```

`FourVectors.jl` is currently unregistered. If you are working with a Julia
version different from the one used by the checked-in manifests in this
repository, install `FourVectors.jl` explicitly before `RemboOnDiet.jl`.

## Docs

- [Documentation website](https://mmikhasenko.github.io/RamboOnDiet.jl/)
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
using FourVectors

masses = [0.93827208816, 0.493677, 0.13957039]
generator = PhaseSpaceGenerator(masses, 2.28646)
point = rand(generator)

weight = phase_space_weight(point)
total = total_momentum(point.momenta)
```
