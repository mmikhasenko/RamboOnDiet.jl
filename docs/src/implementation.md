# Implementation

This page follows the actual code path in `RemboOnDiet.jl` so that the mathematical
description maps directly onto the implementation.

## Public objects

The user-facing types are:

- `PhaseSpaceGenerator`: stores the final-state masses, total four-momentum,
  total invariant mass, threshold, and the implied random-number count,
- `PhaseSpacePoint`: stores the generated four-vectors and the event weight.

Sampling is intentionally idiomatic Julia:

```julia
generator = PhaseSpaceGenerator([0.0, 0.0, 0.0], 2.0)
point = rand(generator)
```

The number of required random inputs is exposed as
`required_random_numbers`, and `generator.random_numbers` is a derived property.

## Data flow in `generate_from_unit_hypercube`

The full deterministic map lives in
[`src/generator.jl`](/Users/mikhailmikhasenko/Documents/JuliaDev.CAT/PhaseSpaceRembo/src/generator.jl).
Given a vector `rs` of length `3n - 4`, the function:

1. checks multiplicity and threshold,
2. computes cumulative tail masses,
3. generates the chain of intermediate cluster masses,
4. performs `n - 1` two-body decays,
5. boosts each decay product back into the lab frame,
6. evaluates the event weight.

The ordering matters because the same `current_cluster` four-vector is both the parent
for the present decay and the boost target for the next one.

## Tail masses and reduced masses

The code first builds

```julia
tail_masses[i] = m_i + m_{i+1} + ... + m_n
```

so that each intermediate cluster mass can be written as

```julia
cluster_masses[i] = reduced_masses[i] + tail_masses[i]
```

This is the massive deformation of the massless chain. In the massless limit,
`tail_masses[i] = 0` and the cluster masses reduce to the reduced masses themselves.

## Solving for `u`

The helper [`solve_mass_parameter`](/Users/mikhailmikhasenko/Documents/JuliaDev.CAT/PhaseSpaceRembo/src/solver.jl)
implements the one-dimensional inverse map from a uniform random number `v` to the
flattening variable `u`. It uses a safeguarded Newton iteration:

- Newton's method provides fast convergence in the interior,
- a bracketing interval `[lo, hi]` keeps the iterate inside `[0, 1]`,
- the fallback is a bisection-like midpoint when a Newton step leaves the bracket.

This solver is cheap, robust, and deterministic, which is useful for testing with a
fixed unit-hypercube input.

## Rotations

The code does not construct explicit rotation matrices. Instead it generates a unit vector

```math
\hat n(\theta,\phi) =
(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)
```

through the helper `unit_direction` in [src/kinematics.jl](/Users/mikhailmikhasenko/Documents/JuliaDev.CAT/PhaseSpaceRembo/src/kinematics.jl). In the cluster rest frame, this is equivalent to taking
a reference momentum on the `z` axis and rotating it by the polar and azimuthal angles.

This direct construction keeps the implementation short and avoids carrying rotation
objects around.

## Boosts

The function `boost` in [src/kinematics.jl](/Users/mikhailmikhasenko/Documents/JuliaDev.CAT/PhaseSpaceRembo/src/kinematics.jl) implements a Lorentz boost defined by the parent four-vector.
If `Q` is the current cluster and `p` is a daughter built in the `Q` rest frame, the code
extracts

```math
\vec\beta = \frac{\vec Q}{Q^0},\qquad
\gamma = \frac{Q^0}{\sqrt{Q^2}}.
```

The boosted momentum is then

```math
\vec p\,' = \vec p +
\left[\frac{\gamma - 1}{\beta^2}(\vec\beta\cdot\vec p) + \gamma p^0\right]\vec\beta,
\qquad
p'^0 = \gamma(p^0 + \vec\beta\cdot\vec p).
```

This is exactly the standard boost of a four-vector by the velocity of the parent cluster.

## Threshold handling

Near threshold, the available kinetic mass goes to zero and the generic recursion becomes
numerically delicate. The package therefore uses a dedicated threshold path:

- in the exact center-of-mass frame, all daughters are returned at rest,
- in a moving frame, each daughter is built at rest and then boosted with the parent.

The threshold weight is returned as zero because the phase-space measure collapses there.

## Why the three-body tests are strong

The package deliberately uses the main sampler for the two-body and three-body validations.
These are the cases where the expected geometry is clearest:

- 2-body: fixed breakup momentum and isotropic angular distribution,
- 3-body massless: flat Dalitz density without weights,
- 3-body massive: flat Dalitz density after weighting by the returned Jacobian.

That makes the three-body Dalitz plots a direct probe of the core mapping rather than of a
special-case implementation.
