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

Sample an aligned phase-space point, redrawing exact-zero mass coordinates.
For two bodies, no random numbers are needed: particle 1 points along +z and
particle 2 along -z in the total rest frame. For three or more bodies, particle 1
points along -z and particle 2 lies in the xz plane with nonnegative x. The weight
includes the integrated overall orientation, as in the ordinary generator.
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

For two bodies, pass an empty vector: particle 1 points along +z and particle 2
along -z in the total rest frame. For `n >= 3`, map `3n - 7` coordinates to an
aligned n-body point.
Above threshold, mass coordinates must be strictly in `(0, 1)` (otherwise a
`DomainError` is thrown); angular coordinates may include both endpoints of `[0, 1]`.
At the total threshold, any finite unit-cube inputs give the threshold configuration.
The first `n - 2` inputs set cluster masses, the next sets `cos(theta_2)`,
and the rest alternate polar and azimuthal coordinates for splittings 3 onward.
Fix `cos(theta_1) = -1`, `phi_1 = 0`, and `phi_2 = 0` in the total rest frame.
For a moving total, boost the completed rest-frame configuration to the lab.
Weights retain the full phase-space normalization with overall rotations integrated
out; use these points for observables invariant under overall rest-frame rotations.
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
