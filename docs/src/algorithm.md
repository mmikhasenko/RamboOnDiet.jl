# Algorithm

`RemboOnDiet.jl` implements the sequential RAMBO-on-diet construction for an `n`-body
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
w = \Phi_n^{(0)}(P^2)\,
\prod_{i=2}^{n-1} \frac{M_i}{\mu_i}
\prod_{i=2}^{n}
\frac{\rho_2(M_{i-1}; m_{i-1}, M_i)}
     {\rho_2(\mu_{i-1}; 0, \mu_i)},
```

with the two-body density

```math
\rho_2(M; m_a, m_b) = \frac{q(M; m_a, m_b)}{4M}.
```

That is why the three-body massless Dalitz plot is flat without weights, while the
massive Dalitz plot becomes flat only after filling with `weights = phase_space_weight(point)`.
