struct PhaseSpaceGenerator{T,N,M<:AbstractVector{T}}
    masses::M
    total::FourVector{T}
    invariant_mass::T
    threshold::T
end

Base.length(::PhaseSpaceGenerator{T,N}) where {T,N} = N

function Base.getproperty(generator::PhaseSpaceGenerator{T,N}, name::Symbol) where {T,N}
    name === :random_numbers && return required_random_numbers(N)
    return getfield(generator, name)
end

function Base.propertynames(::PhaseSpaceGenerator, private::Bool = false)
    names = (:masses, :total, :invariant_mass, :threshold, :random_numbers)
    return private ? names : names
end

function PhaseSpaceGenerator(
    masses::AbstractVector{<:Real},
    total::FourVector{T};
    dalitz_grid::Int = 4097,
) where {T}
    massesT = T.(masses)
    N = length(massesT)
    invariant_mass = mass(total)
    threshold = sum(massesT)
    return PhaseSpaceGenerator{T,N,typeof(massesT)}(massesT, total, invariant_mass, threshold)
end

function PhaseSpaceGenerator(
    masses::NTuple{N,<:Real},
    total::FourVector{T};
    dalitz_grid::Int = 4097,
) where {N,T}
    return PhaseSpaceGenerator(SVector{N,T}(masses), total; dalitz_grid)
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

function PhaseSpaceGenerator(
    masses::NTuple{N,<:Real},
    sqrt_s::Real;
    dalitz_grid::Int = 4097,
) where {N}
    T = float(promote_type(eltype(masses), typeof(sqrt_s)))
    total = FourVector(zero(T), zero(T), zero(T); E = T(sqrt_s))
    return PhaseSpaceGenerator(SVector{N,T}(masses), total; dalitz_grid)
end
