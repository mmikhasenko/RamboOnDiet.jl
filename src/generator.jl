required_random_numbers(n::Integer) = 3 * n - 4

# Draw the whole cube first; only mass coordinates need zero resampling.
function resample_mass_zeros!(rng::AbstractRNG, rs::AbstractVector, n_mass::Integer)
    for i = 1:n_mass
        while iszero(rs[i])
            rs[i] = rand(rng, eltype(rs))
        end
    end
    return rs
end

"""
    generate_point([rng], generator::PhaseSpaceGenerator)
    generate_point([rng], masses, sqrt_s)

Draw one `n`-body phase-space event. Masses and `sqrt_s` use one energy unit.
`rand([rng], generator)` calls this function.

Build the parent with [`PhaseSpaceGenerator`](@ref), or pass masses and
`sqrt_s` for a parent at rest. The result is a [`PhaseSpacePoint`](@ref).
Read the daughters from `point.momenta` or `point[i]`, then the Monte Carlo
weight with [`phase_space_weight`](@ref).

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
point = generate_point(generator)
point.momenta
```

## Related

- [`PhaseSpaceGenerator`](@ref) — masses and the total four-momentum
- [`PhaseSpacePoint`](@ref) — momenta and weight returned by this draw
- [`phase_space_weight`](@ref) — weight to accumulate for this event
- [`generate_from_unit_hypercube`](@ref) — run the same map at fixed coordinates
"""
function generate_point(rng::AbstractRNG, generator::PhaseSpaceGenerator{T}) where {T}
    n = length(generator)
    rs = rand(rng, T, required_random_numbers(n))
    resample_mass_zeros!(rng, rs, n - 2)
    return generate_from_unit_hypercube(rs, generator)
end

generate_point(generator::PhaseSpaceGenerator) =
    generate_point(Random.default_rng(), generator)
Base.rand(rng::AbstractRNG, generator::PhaseSpaceGenerator) = generate_point(rng, generator)
Base.rand(generator::PhaseSpaceGenerator) = generate_point(generator)

function generate_momenta(
    rng::AbstractRNG,
    masses::AbstractVector{<:Real},
    total::FourVector{T},
) where {T}
    generator = PhaseSpaceGenerator(masses, total)
    return generate_point(rng, generator)
end

function generate_point(rng::AbstractRNG, masses::AbstractVector{<:Real}, sqrt_s::Real)
    generator = PhaseSpaceGenerator(masses, sqrt_s)
    return generate_point(rng, generator)
end

generate_point(masses::AbstractVector{<:Real}, sqrt_s::Real) =
    generate_point(Random.default_rng(), masses, sqrt_s)
generate_momenta(masses::AbstractVector{<:Real}, total::FourVector) =
    generate_momenta(Random.default_rng(), masses, total)

"""
    generate_from_unit_hypercube(rs, generator::PhaseSpaceGenerator)
    generate_from_unit_hypercube(rs, masses, total::FourVector)

Map one unit-cube vector to one phase-space event. Masses and the total
energy use one energy unit.

Build the parent with [`PhaseSpaceGenerator`](@ref). `rs` has length
`generator.random_numbers` (`3n - 4`): the first `n - 2` entries are mass
coordinates and lie in `(0, 1)` above threshold, and the rest are angles in
`[0, 1]`. The map returns a [`PhaseSpacePoint`](@ref). Read `point.weight`,
or call [`phase_space_weight`](@ref), before using the daughters.

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
point = generate_from_unit_hypercube([0.3, 0.2, 0.7, 0.4, 0.9], generator)
point.weight
```

## Related

- [`PhaseSpaceGenerator`](@ref) — masses, total momentum, and the coordinate count
- [`generate_point`](@ref) — draw those coordinates at random
- [`generate_aligned_from_unit_hypercube`](@ref) — pass the `3n - 7` aligned coordinates
- [`PhaseSpacePoint`](@ref) — momenta and weight for this coordinate
- [`phase_space_weight`](@ref) — Monte Carlo weight of this coordinate
"""
function generate_from_unit_hypercube(
    rs::AbstractVector{<:Real},
    generator::PhaseSpaceGenerator{T},
) where {T}
    return generate_from_unit_hypercube(rs, generator.masses, generator.total)
end

function generate_from_unit_hypercube(
    rs::AbstractVector{<:Real},
    masses::AbstractVector{<:Real},
    total::FourVector{T},
) where {T}
    n = length(masses)
    n >= 2 || throw(ArgumentError("At least two outgoing particles are required."))
    length(rs) == required_random_numbers(n) || throw(
        ArgumentError(
            "Expected $(required_random_numbers(n)) random numbers for $n particles, got $(length(rs)).",
        ),
    )

    all(r -> isfinite(r) && 0 <= r <= 1, rs) ||
        throw(ArgumentError("Coordinates must be finite and in [0, 1]."))

    massesT = T.(masses)
    Mtot = LorentzVectorBase.mass(total)
    total_mass = sum(massesT)
    Mtot + sqrt(eps(T)) < total_mass &&
        throw(ArgumentError("Total invariant mass $Mtot is below threshold $total_mass."))

    if isapprox(Mtot, total_mass; atol = sqrt(eps(T)), rtol = sqrt(eps(T)))
        return threshold_configuration(massesT, total)
    end

    # Above threshold, the mass map is defined only on the open interval.
    # Check in the working precision too: conversion can round an input to 0 or 1.
    for i = 1:(n-2)
        r = T(rs[i])
        0 < r < 1 || throw(
            DomainError(rs[i], "Mass coordinate $i must lie strictly between 0 and 1."),
        )
    end

    tail_masses = similar(massesT)
    acc = zero(T)
    for i = n:-1:1
        acc += massesT[i]
        tail_masses[i] = acc
    end

    cluster_masses = Vector{T}(undef, n)
    reduced_masses = Vector{T}(undef, n)
    cluster_masses[1] = Mtot
    reduced_masses[1] = Mtot - total_mass
    reduced_masses[1] < zero(T) &&
        throw(ArgumentError("Reduced mass is negative; kinematics are invalid."))

    idx = 1
    for i = 2:(n-1)
        power = n - i
        u = solve_mass_parameter(T(rs[idx]), power)
        idx += 1
        reduced_masses[i] = reduced_masses[i-1] * sqrt(u)
        cluster_masses[i] = reduced_masses[i] + tail_masses[i]
    end
    cluster_masses[n] = massesT[n]
    reduced_masses[n] = zero(T)

    momenta = Vector{FourVector{T}}(undef, n)
    current_cluster = total

    for i = 1:(n-1)
        cosθ = T(2 * rs[idx] - 1)
        ϕ = T(2π * rs[idx+1])
        idx += 2

        child_mass = cluster_masses[i+1]
        q = breakup_momentum(cluster_masses[i], massesT[i], child_mass)
        direction = unit_direction(cosθ, ϕ)
        daughter_cm = fourvector(q * direction, sqrt(q^2 + massesT[i]^2))
        child_cm = fourvector(-q * direction, sqrt(q^2 + child_mass^2))

        momenta[i] = boost(daughter_cm, current_cluster)
        current_cluster = boost(child_cm, current_cluster)
    end
    momenta[n] = current_cluster

    # The auxiliary massless phase space starts at the available kinetic mass Δ.
    base_weight = T(massless_phase_space_volume(reduced_masses[1]^2, n))
    if all(iszero, massesT)
        return PhaseSpacePoint(momenta, base_weight)
    end

    jacobian = one(T)
    for i = 2:(n-1)
        jacobian *= cluster_masses[i] / reduced_masses[i]
    end
    for i = 2:n
        jacobian *= two_body_density(cluster_masses[i-1], massesT[i-1], cluster_masses[i])
        jacobian /= two_body_density(reduced_masses[i-1], zero(T), reduced_masses[i])
    end

    return PhaseSpacePoint(momenta, base_weight * jacobian)
end

function threshold_configuration(masses::AbstractVector{T}, total::FourVector{T}) where {T}
    Mtot = LorentzVectorBase.mass(total)
    if isapprox(
        LorentzVectorBase.spatial_magnitude(total),
        zero(T);
        atol = sqrt(eps(T)),
        rtol = sqrt(eps(T)),
    )
        momenta = [FourVector(zero(T), zero(T), zero(T); E = m) for m in masses]
        return PhaseSpacePoint(momenta, zero(T))
    end
    γ = LorentzVectorBase.E(total) / Mtot
    momenta = Vector{FourVector{T}}(undef, length(masses))
    for (i, m) in enumerate(masses)
        rest_p = FourVector(zero(T), zero(T), zero(T); E = m)
        momenta[i] = boost(rest_p, total)
        # Guard against accumulated roundoff if the user lands exactly on threshold.
        momenta[i] = fourvector(spatial(momenta[i]), γ * m)
    end
    return PhaseSpacePoint(momenta, zero(T))
end

# Internal coordinate count for aligned phase space.
function required_aligned_random_numbers(n::Integer)
    n >= 2 || throw(ArgumentError("Aligned generation requires at least two particles."))
    n == 2 && return 0
    return 3 * n - 7
end

"""
    generate_aligned_point([rng], generator)
    generate_aligned_point([rng], masses, total_or_sqrt_s)

Draw one event with the overall rotation fixed. Masses and the total energy
use one energy unit.

The parent is a [`PhaseSpaceGenerator`](@ref). Two bodies draw no random
numbers: in the total rest frame, particle 1 points along +z and particle 2
along -z. For three or more bodies, particle 1 points along -z and particle 2
lies in the xz plane with nonnegative x. The sampler draws `3n - 7`
coordinates and redraws an exact-zero mass coordinate. The returned
[`PhaseSpacePoint`](@ref) uses the ordinary normalization, with that rotation
integrated into [`phase_space_weight`](@ref).

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
point = generate_aligned_point(generator)
point[1]
```

## Related

- [`PhaseSpaceGenerator`](@ref) — masses and the total four-momentum
- [`generate_aligned_from_unit_hypercube`](@ref) — supply the reduced coordinates yourself
- [`generate_point`](@ref) — sample the overall rotation as well
- [`PhaseSpacePoint`](@ref) — aligned momenta and the stored weight
- [`phase_space_weight`](@ref) — orientation-integrated weight of this draw
"""
function generate_aligned_point(
    rng::AbstractRNG,
    generator::PhaseSpaceGenerator{T},
) where {T}
    n = length(generator)
    rs = rand(rng, T, required_aligned_random_numbers(n))
    resample_mass_zeros!(rng, rs, n - 2)
    return generate_aligned_from_unit_hypercube(rs, generator)
end

generate_aligned_point(generator::PhaseSpaceGenerator) =
    generate_aligned_point(Random.default_rng(), generator)
generate_aligned_point(rng::AbstractRNG, masses, total_or_sqrt_s) =
    generate_aligned_point(rng, PhaseSpaceGenerator(masses, total_or_sqrt_s))
generate_aligned_point(masses, total_or_sqrt_s) =
    generate_aligned_point(Random.default_rng(), masses, total_or_sqrt_s)

"""
    generate_aligned_from_unit_hypercube(rs, generator)
    generate_aligned_from_unit_hypercube(rs, masses, total_or_sqrt_s)

Map reduced unit-cube coordinates to one aligned event. Masses and the total
energy use one energy unit.

Build a [`PhaseSpaceGenerator`](@ref) first. For two bodies, pass an empty
vector: particle 1 points along +z and particle 2 along -z in the total rest
frame. For `n ≥ 3`, `rs` has length `3n - 7`. The first `n - 2` entries set
the cluster masses and lie in `(0, 1)` above threshold. The next entry is the
polar coordinate of particle 2, and the rest are polar/azimuthal pairs for the
later splittings. Angles may be 0 or 1. A moving `generator.total` is applied
after the rest-frame configuration is built. The result is a [`PhaseSpacePoint`](@ref).

```julia
using RamboOnDiet
two_body = PhaseSpaceGenerator([0.3, 0.5], 2.0)
generate_aligned_from_unit_hypercube(Float64[], two_body)
```

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
point = generate_aligned_from_unit_hypercube([0.4, 0.6], generator)
point.momenta
```

## Related

- [`PhaseSpaceGenerator`](@ref) — masses and the total four-momentum being decayed
- [`generate_aligned_point`](@ref) — draw the reduced coordinates at random
- [`generate_from_unit_hypercube`](@ref) — keep the overall rotation and pass `3n - 4` coordinates
- [`PhaseSpacePoint`](@ref) — momenta and weight for these coordinates
- [`phase_space_weight`](@ref) — weight with the overall rotation integrated out
"""
function generate_aligned_from_unit_hypercube(
    rs::AbstractVector{<:Real},
    generator::PhaseSpaceGenerator{T},
) where {T}
    n = length(generator)
    expected = required_aligned_random_numbers(n)
    length(rs) == expected || throw(
        ArgumentError(
            "Expected $expected coordinates for $n aligned particles, got $(length(rs)).",
        ),
    )
    all(r -> isfinite(r) && 0 <= r <= 1, rs) ||
        throw(ArgumentError("Aligned coordinates must be finite and in [0, 1]."))

    n == 2 && return generate_from_unit_hypercube(T[1, 0], generator)

    # Restore the three fixed angular inputs in the ordinary cube ordering.
    full = Vector{T}(undef, required_random_numbers(n))
    full[1:(n-2)] .= rs[1:(n-2)]
    full[n-1] = zero(T) # cos(theta_1) = -1
    full[n] = zero(T)   # phi_1 = 0
    full[n+1] = rs[n-1] # free polar coordinate for particle 2
    full[n+2] = zero(T) # phi_2 = 0
    full[(n+3):end] .= rs[n:end]

    rest_total = FourVector(zero(T), zero(T), zero(T); E = generator.invariant_mass)
    point = generate_from_unit_hypercube(full, generator.masses, rest_total)
    momenta = [boost(p, generator.total) for p in point.momenta]
    return PhaseSpacePoint(momenta, point.weight)
end

generate_aligned_from_unit_hypercube(rs::AbstractVector{<:Real}, masses, total_or_sqrt_s) =
    generate_aligned_from_unit_hypercube(rs, PhaseSpaceGenerator(masses, total_or_sqrt_s))
