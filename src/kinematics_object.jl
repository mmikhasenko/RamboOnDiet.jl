struct PhaseSpaceGenerator{T}
    masses::Vector{T}
    total::FourVector{T}
    invariant_mass::T
    threshold::T
    random_numbers::Int
end

function PhaseSpaceGenerator(
    masses::AbstractVector{<:Real},
    total::FourVector{T};
    dalitz_grid::Int = 4097,
) where {T}
    massesT = T.(masses)
    invariant_mass = mass(total)
    threshold = sum(massesT)
    return PhaseSpaceGenerator{T}(massesT, total, invariant_mass, threshold, required_random_numbers(length(massesT)))
end

function PhaseSpaceGenerator(
    masses::AbstractVector{<:Real},
    sqrt_s::Real;
    dalitz_grid::Int = 4097,
)
    T = float(promote_type(eltype(masses), typeof(sqrt_s)))
    total = FourVector(zero(T), zero(T), zero(T); E = T(sqrt_s))
    return PhaseSpaceGenerator(masses, total; dalitz_grid)
end
