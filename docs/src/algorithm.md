# Algorithm

`RamboOnDiet.jl` implements the sequential RAMBO-on-diet construction for an `n`-body
phase-space point with total four-momentum `P` and outgoing masses `m_1, \ldots, m_n`.
The public entry point is `PhaseSpaceGenerator`, but the core map is easiest to
read as a sequence of deterministic steps fed by `3n - 4` uniform random numbers.

## Goal

We want to sample the Lorentz-invariant phase-space measure

```math
d\Phi_n(P; p_1,\ldots,p_n) =
\left[\prod_{i=1}^n d^4 p_i \, \delta^+(p_i^2 - m_i^2)\right]
\delta^{(4)}\!\left(P - \sum_{i=1}^n p_i\right).
```

The RAMBO-on-diet strategy factorizes this into a chain of two-body decays,

```math
P = Q_1 \to p_1 + Q_2,\quad
Q_2 \to p_2 + Q_3,\quad \ldots,\quad
Q_{n-1} \to p_{n-1} + p_n,
```

where the intermediate cluster masses `M_i = \sqrt{Q_i^2}` are generated from the
unit hypercube.

This is the same general strategy introduced for massless phase space by Kleiss, Stirling,
and Ellis in the original RAMBO paper and then recast into the minimal `3n - 4`
parameterization by Plätzer's RAMBO-on-diet construction. In this package,
`solve_mass_parameter` numerically inverts the one-dimensional cumulative map that appears
in that RAMBO-on-diet mass parametrization.

References:

- R. Kleiss, W. J. Stirling, and S. D. Ellis, "A new Monte Carlo treatment of multiparticle
  phase space at high energies", Comput. Phys. Commun. 40 (1986) 359-373,
  [doi:10.1016/0010-4655(86)90119-0](https://doi.org/10.1016/0010-4655(86)90119-0).
- S. Plätzer, "RAMBO on diet", [arXiv:1308.2922](https://arxiv.org/abs/1308.2922).

## Random variables

For `n` particles the generator consumes

```math
3n - 4 = (n - 2) + 2(n - 1)
```

uniform random numbers:

- `n - 2` numbers determine the chain of intermediate masses,
- `2(n - 1)` numbers determine the decay angles of the sequential two-body decays.

The last intermediate mass is fixed by kinematics, so only `n - 2` mass-like variables
are needed.

## Mapping of the `u` variables

Let

```math
\Delta = \sqrt{P^2} - \sum_{i=1}^n m_i
```

be the kinetic mass available above threshold. The implementation stores the cumulative
"tail" masses

```math
\tau_i = \sum_{k=i}^n m_k
```

and defines reduced masses

```math
\mu_1 = \Delta,\qquad
M_i = \mu_i + \tau_i.
```

For each `i = 2,\ldots,n-1`, a uniform random number `r_i` is converted into a parameter
`u_i \in [0,1]` by solving

```math
r_i = (a_i + 1) u_i^{a_i} - a_i u_i^{a_i + 1},
\qquad a_i = n - i.
```

This is the flattened massless RAMBO-on-diet distribution for the sequential mass ratios.
The reduced masses are then propagated as

```math
\mu_i = \mu_{i-1} \sqrt{u_i},
\qquad
M_i = \mu_i + \tau_i.
```

The square root is important: the recursion is naturally written in terms of invariant
masses squared, so the masses themselves inherit a `\sqrt{u_i}` factor.

## Angles and two-body decays

For each decay step `Q_i \to p_i + Q_{i+1}`, the generator takes two more uniform random
numbers and maps them to

```math
\cos\theta_i = 2r - 1,\qquad \phi_i = 2\pi r'.
```

In the rest frame of the decaying cluster `Q_i`, the daughter momentum magnitude is

```math
q_i = \frac{\sqrt{\lambda(M_i^2, m_i^2, M_{i+1}^2)}}{2 M_i},
```

with the Källén function

```math
\lambda(x,y,z) = x^2 + y^2 + z^2 - 2(xy + yz + zx).
```

The spatial direction is built from `(\cos\theta_i,\phi_i)`, the visible daughter is
placed at `+\vec q_i`, and the next cluster at `-\vec q_i`.

## Boost chain

The sequential decays are generated in their own cluster rest frames, then each daughter
and child cluster is boosted into the current lab frame. Starting from `Q_1 = P`,

1. construct `p_i` and `Q_{i+1}` in the `Q_i` rest frame,
2. boost both four-vectors by the velocity of `Q_i`,
3. keep the boosted `Q_{i+1}` as the parent of the next decay.

This guarantees exact four-momentum conservation at every step, up to floating-point
roundoff.

## Event weight

The returned `PhaseSpacePoint` stores a weight:

- for `m_i = 0`, the weight is constant and equal to the massless phase-space volume,
- for `m_i \neq 0`, the same massless map is corrected by the RAMBO-on-diet Jacobian.

In the implementation this is written as

```math
w = \Phi_n^{(0)}(\Delta^2)\,
\prod_{i=2}^{n-1} \frac{M_i}{\mu_i}
\prod_{i=2}^{n}
\frac{\rho_2(M_{i-1}; m_{i-1}, M_i)}
     {\rho_2(\mu_{i-1}; 0, \mu_i)},
```

The prefactor is the massless volume at the reduced mass `Δ`, because the auxiliary
mass chain starts at `μ_1 = Δ`. For `n > 2`, it equals the volume at `P²` only when
all daughter masses vanish. Using `P²` for massive daughters would overestimate every weight by
`(√(P²)/Δ)^(2n - 4)`.

The two-body density is

```math
\rho_2(M; m_a, m_b) = \frac{q(M; m_a, m_b)}{4M}.
```

In `src/generator.jl`, this factorization appears as two explicit products: one loop over `i = 2:(n - 1)`
for the `M_i / \mu_i` terms and one loop over `i = 2:n` for the ratio of massive to
massless two-body densities.

That is why the three-body massless Dalitz plot is flat without weights, while the
massive Dalitz plot becomes flat only after filling with `weights = phase_space_weight(point)`.

## Coordinate ordering: masses first, then angles

For **n final particles there are n - 1 two-body splittings**. The `3n - 4`
unit-cube inputs are ordered as follows:

- The **first n - 2 coordinates** determine the intermediate cluster masses.
- The **last 2(n - 1) coordinates** alternate polar and azimuthal coordinates,
  one pair per splitting, in that splitting's parent rest frame.

These inputs map to masses and angles; they are not themselves dimensional masses
or angles. Polar inputs are uniform in `cos(theta)`, not in `theta`.
The mass map is nonlinear and coupled: each cluster mass sets the range for the next.

In the tables below, `u_1, u_2, ...` name the **unit-cube inputs** (called `rs` in
code), not the derived mass-ratio variables `u_i` in the earlier mass recursion.
Let `Delta = M - sum(masses)` with `M = sqrt(P²)`.

### Three bodies

The chain is `P -> p1 + Q23`, then `Q23 -> p2 + p3`.

| Input | Meaning | Map |
|:--|:--|:--|
| `u_1` | Pair mass `m23` | `m23 = m2 + m3 + Delta * sqrt(1 - sqrt(1 - u_1))` |
| `u_2` | Polar angle of `p1` in the total rest frame | `cos(theta_1) = 2u_2 - 1` |
| `u_3` | Azimuth of `p1` in the total rest frame | `phi_1 = 2pi * u_3` |
| `u_4` | Polar angle of `p2` in the `Q23` rest frame | `cos(theta_2) = 2u_4 - 1` |
| `u_5` | Azimuth of `p2` in the `Q23` rest frame | `phi_2 = 2pi * u_5` |

Thus there is one mass coordinate and four angular coordinates. The other independent
Dalitz invariant is determined by the relative decay directions.

### Four bodies

The chain is `P -> p1 + Q234`, then `Q234 -> p2 + Q34`, then `Q34 -> p3 + p4`.
Let `v` be the unique root in `[0, 1]` of `3v² - 2v³ = u_1`.

| Input | Meaning | Map |
|:--|:--|:--|
| `u_1` | Cluster mass `m234` | `m234 = m2 + m3 + m4 + Delta * sqrt(v)` |
| `u_2` | Pair mass `m34` | `m34 = m3 + m4 + Delta * sqrt(v) * sqrt(1 - sqrt(1 - u_2))` |
| `u_3` | Polar angle of `p1` in the total rest frame | `cos(theta_1) = 2u_3 - 1` |
| `u_4` | Azimuth of `p1` in the total rest frame | `phi_1 = 2pi * u_4` |
| `u_5` | Polar angle of `p2` in the `Q234` rest frame | `cos(theta_2) = 2u_5 - 1` |
| `u_6` | Azimuth of `p2` in the `Q234` rest frame | `phi_2 = 2pi * u_6` |
| `u_7` | Polar angle of `p3` in the `Q34` rest frame | `cos(theta_3) = 2u_7 - 1` |
| `u_8` | Azimuth of `p3` in the `Q34` rest frame | `phi_3 = 2pi * u_8` |

The axes are those used by the rotation-free boosts in the implementation. The
polar angles are not automatically helicity angles along each cluster's flight direction.

## Aligned four-vectors

For `n >= 3`, overall spatial orientation can be fixed by choosing:

- `p1` along **-z** in the total rest frame: `cos(theta_1) = -1` and `phi_1 = 0`;
- `phi_2 = 0` in the `Q_{2...n}` rest frame, placing `p2` in the **xz plane with x >= 0**.

The conventional polar angle for -z is `pi`; using `-pi` gives the same direction.
These are three fixed angular coordinates, but they are **not three consecutive
entries**: `theta_2` remains free. Since the first cluster travels along +z,
its boost also leaves `p2` in the positive-x half of the xz plane in the total rest frame.

The reduced input belongs to `[0, 1]^(3n - 7)`, with ordering:

| Reduced inputs | Meaning |
|:--|:--|
| First `n - 2` | The same intermediate mass coordinates as the ordinary map |
| Next one | Polar coordinate for `p2`: `cos(theta_2) = 2u - 1` |
| Remaining `2(n - 3)` | Polar/azimuthal pairs for splittings `3, ..., n - 1` |

For three bodies this is `[mass23, polar2]`; for four bodies it is
`[mass234, mass34, polar2, polar3, azimuth3]`, where each entry is a unit-cube
coordinate mapped as above.

```@example aligned
using Random, RamboOnDiet

generator = PhaseSpaceGenerator([0.1, 0.2, 0.3, 0.4], 3.0)
required_aligned_random_numbers(4) # 5

# Random sampling, with an optional explicit RNG.
point = generate_aligned_point(MersenneTwister(7), generator)

# Deterministic map from the reduced unit cube.
point = generate_aligned_from_unit_hypercube([0.3, 0.6, 0.4, 0.7, 0.2], generator)
phase_space_weight(point)
```

Both methods also accept `(masses, sqrt_s)` or `(masses, total_fourvector)` in
place of the generator. For a moving total, the entire aligned configuration is
first built in the total rest frame and then boosted to the lab. The alignment
conditions therefore apply in the total rest frame, not generally in the lab.
At threshold or collinear boundaries, some directions are undefined and the
orientation convention is degenerate.

**Weights include the integrated overall orientation.** They have the same
normalization as ordinary points; do not multiply them by another rotation-volume
factor. Aligned samples are appropriate for observables invariant under overall
rotations in the total rest frame. To model orientation-dependent observables,
use the ordinary sampler or restore a uniformly sampled overall rotation.

```@docs
required_aligned_random_numbers
generate_aligned_point
generate_aligned_from_unit_hypercube
```
