required_random_numbers(n::Integer) = 3 * n - 4

function generate_point(rng::AbstractRNG, generator::PhaseSpaceGenerator{T,N}) where {T,N}
    rs = SVector{required_random_numbers(N),T}(
        ntuple(_ -> rand(rng, T), Val(required_random_numbers(N))),
    )
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

    massesT = T.(masses)
    Mtot = LorentzVectorBase.mass(total)
    total_mass = sum(massesT)
    Mtot + sqrt(eps(T)) < total_mass &&
        throw(ArgumentError("Total invariant mass $Mtot is below threshold $total_mass."))

    if isapprox(Mtot, total_mass; atol = sqrt(eps(T)), rtol = sqrt(eps(T)))
        return threshold_configuration(massesT, total)
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

    base_weight = T(massless_phase_space_volume(Mtot^2, n))
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
