```@meta
EditURL = "../../four-body-mass-window.jl"
```

# Four-body sampling in a three-pion mass window

A 190 GeV/c π⁻ beam on a proton at rest produces `[recoil, π⁻, π⁺, π⁻]`.
Particle 1 is the recoil. In the
[four-body coordinate table](../algorithm.md#Four-bodies),
`u₁` is the mass of the cluster formed by particles 2, 3, and 4, so here it
is the three-pion mass. A different mass order windows a different cluster.

Drawing `u₁` from the image of `[1.02, 1.7] GeV` keeps every event inside the
window. The weight is multiplied by the width of that coordinate interval.

The weighted sample estimates the restricted measure

```math
d\Phi_4
= \delta^{4}\!\left(P - \sum_i p_i\right)
\prod_i \frac{d^{3}p_i}{2 E_i},
```

without factors of `2π`, the same convention as `phase_space_weight`.
Absolute massive weights in this normalization require RamboOnDiet v0.5.0 or later.

## Beam and window

The total energy is fixed by the beam momentum. For four bodies the first cube
input is the three-pion mass through

```math
3v^2 - 2v^3 = u_1,
\qquad
m_{3\pi} = 3m_\pi + \Delta\sqrt{v},
\qquad
\Delta = W - m_p - 3m_\pi.
```

Inverting that map gives the coordinate interval the sampler is allowed to draw.

````julia
using RamboOnDiet
using FourVectors
using Random
using LinearAlgebra
using Printf

const mpi = 0.13957039
const mp = 0.93827208816
const beam_momentum = 190.0
const W = sqrt(mp^2 + mpi^2 + 2mp * sqrt(beam_momentum^2 + mpi^2))
const mass_min, mass_max = 1.02, 1.7
const kinetic_mass = W - mp - 3mpi

function mass_coordinate(m)
    z = ((m - 3mpi) / kinetic_mass)^2
    return z^2 * (3 - 2z)
end

const coordinate_min = mass_coordinate(mass_min)
const coordinate_max = mass_coordinate(mass_max)
const coordinate_width = coordinate_max - coordinate_min
const generator = PhaseSpaceGenerator([mp, mpi, mpi, mpi], W)
````

The collision energy is **18.906172 GeV**. The mass window maps onto `u₁ ∈ [4.131e-06, 8.494e-05]`. That interval is **0.008081%** of the cube, so rejection sampling would keep one event in **12374**.

## One event inside the window

`generate_from_unit_hypercube` takes each coordinate as uniform on `(0, 1)`,
where the density in `u₁` is 1. A uniform draw on the subinterval `[u_a, u_b]`
has density `1/(u_b - u_a)` instead. Multiplying `phase_space_weight` by
`u_b - u_a` cancels that density and restores the phase-space measure on the
window. The events stay weighted: the factor is the Jacobian of the
restriction, not an accept/reject decision.

The second coordinate is the remaining mass input and has to stay in the open
unit interval. The daughters come out as `[recoil, π⁻, π⁺, π⁻]`.

````julia
function restricted_event(rng)
    coordinates = rand(rng, 8)
    coordinates[1] = coordinate_min + coordinate_width * coordinates[1]
    # Mass coordinates are defined on the open unit interval.
    while iszero(coordinates[2])
        coordinates[2] = rand(rng)
    end
    event = generate_from_unit_hypercube(coordinates, generator)
    weight = coordinate_width * phase_space_weight(event)
    return event, weight, coordinates[1]
end
````

## Reference spectrum

The same restricted volume has an independent factorization. At three-pion
mass `m` the four-body measure splits into the recoil and the three-pion decay,

```math
\frac{d\Phi_4}{dm} = 2m \, \Phi_2(W; m_p, m) \, \Phi_3(m; m_\pi, m_\pi, m_\pi).
```

The three-body volume is an integral of two two-body volumes over the pair
mass squared. The two-body volume below is angle-integrated and uses the
generator's convention, without factors of `2π`.
Gauss–Legendre rules of 64 and 128 nodes check that this integral has settled.

````julia
phi2(M, a, b) =
    π * sqrt(max(0.0, (M - a - b) * (M + a + b) * (M - a + b) * (M + a - b))) / (2M^2)

function gauss_legendre(n)
    off_diagonal = [k / sqrt(4k^2 - 1) for k = 1:(n-1)]
    decomposition = eigen(SymTridiagonal(zeros(n), off_diagonal))
    return decomposition.values, 2 .* decomposition.vectors[1, :] .^ 2
end

function integrate(f, lower, upper, rule)
    nodes, weights = rule
    midpoint = (lower + upper) / 2
    half_width = (upper - lower) / 2
    return half_width *
           sum(weights[k] * f(midpoint + half_width * nodes[k]) for k in eachindex(nodes))
end

function spectrum(m, rule)
    phi3 = integrate((2mpi)^2, (m - mpi)^2, rule) do t
        phi2(m, mpi, sqrt(t)) * phi2(sqrt(t), mpi, mpi)
    end
    return 2m * phi2(W, mp, m) * phi3
end

function reference_bins(edges, rule)
    return [
        integrate(m -> spectrum(m, rule), edges[b], edges[b+1], rule) for
        b = 1:(length(edges)-1)
    ]
end

````

## Sample

One hundred thousand events, seed `20261021`. The weights estimate the
absolute integral of $d\Phi_4$ in each mass bin.

````julia
function sample_window(; n = 100_000, seed = 20261021)
    rng = MersenneTwister(seed)
    edges = collect(range(mass_min, mass_max; length = 11))
    counts, sums, sums2 = zeros(Int, 10), zeros(10), zeros(10)
    minmass, maxmass = Inf, -Inf
    shell_error = conservation_error = map_error = angle_error = 0.0
    positive_weights = finite_angles = true
    for _ = 1:n
        event, weight, u1 = restricted_event(rng)
        momenta = event.momenta
        parent = sum(momenta)
        m3pi = mass(momenta[2] + momenta[3] + momenta[4])
        discrepancy, finite = angle_discrepancy(momenta, m3pi)
        angle_error = max(angle_error, discrepancy)
        finite_angles &= finite
        positive_weights &= isfinite(weight) && weight > 0
        minmass, maxmass = min(minmass, m3pi), max(maxmass, m3pi)
        map_error = max(map_error, abs(mass_coordinate(m3pi) - u1))
        conservation_error = max(
            conservation_error,
            abs(parent.E - W),
            abs(parent.px),
            abs(parent.py),
            abs(parent.pz),
        )
        shell_error = max(
            shell_error,
            maximum(abs(mass(p) - mi) for (p, mi) in zip(momenta, (mp, mpi, mpi, mpi))),
        )
        bin = clamp(searchsortedlast(edges, m3pi), 1, 10)
        counts[bin] += 1
        sums[bin] += weight
        sums2[bin] += weight^2
    end
    return (;
        n,
        seed,
        edges,
        counts,
        sums,
        sums2,
        minmass,
        maxmass,
        shell_error,
        conservation_error,
        map_error,
        angle_error,
        positive_weights,
        finite_angles,
    )
end
````

## Result

The table compares the weighted Monte Carlo integral in each bin with the
128-point quadrature. Pulls use the Monte Carlo standard error. A second rule,
with 64 points, checks the quadrature itself.

Raw counts rise toward the upper edge of the window because the map from `u₁`
to the three-pion mass is nonlinear. Those counts follow the restricted cube
map. The weighted fraction is the column that should match the reference
fraction.

````julia
sample = sample_window()
coarse = reference_bins(sample.edges, gauss_legendre(64))
exact = reference_bins(sample.edges, gauss_legendre(128))
````

Seed `20261021`, **100,000** events. Every three-pion mass lies in
[1.020036920, 1.699993129] GeV.

The restricted volume is **5.1253743 ± 0.00422**.
Quadrature gives **5.1311734** (pull **-1.37**).
The effective sample size is **93,647** events (93.65%).
The 64-point rule matches the 128-point rule to a relative **3.91e-06**.

Momentum is conserved to **2.37e-12 GeV**, masses sit on shell to **2.30e-13 GeV**,
and the cube map inverts to **7.10e-15**. Decay angles agree with the Cartesian boosts to **6.36e-11 rad**.
Unweighted counts follow the cube map; the largest of those pulls is **1.61**.

| m₃π (GeV) | counts | weighted fraction | reference fraction | pull |
|:--|--:|--:|--:|--:|
| 1.020–1.088 | 2722 | 0.035521 | 0.035493 | -0.02 |
| 1.088–1.156 | 3637 | 0.044501 | 0.045297 | -1.13 |
| 1.156–1.224 | 4874 | 0.056680 | 0.056634 | -0.02 |
| 1.224–1.292 | 6372 | 0.070078 | 0.069600 | +0.46 |
| 1.292–1.360 | 7992 | 0.085021 | 0.084290 | +0.67 |
| 1.360–1.428 | 9819 | 0.100368 | 0.100797 | -0.55 |
| 1.428–1.496 | 11964 | 0.118101 | 0.119216 | -1.20 |
| 1.496–1.564 | 14438 | 0.138514 | 0.139640 | -1.17 |
| 1.564–1.632 | 17408 | 0.162246 | 0.162161 | -0.09 |
| 1.632–1.700 | 20774 | 0.188970 | 0.186871 | +1.57 |

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

