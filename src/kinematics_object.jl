"""
    PhaseSpaceGenerator(masses, sqrt_s)
    PhaseSpaceGenerator(masses, total::FourVector)

Parent process for an `n`-body RAMBO-on-diet sampler.

`masses` are the outgoing masses and `sqrt_s` is the center-of-mass energy,
in one common energy unit. Pass a `FourVector` as `total` when the parent is
moving; its components use that same unit. The object stores `generator.masses`,
`generator.total`, `generator.invariant_mass`, and `generator.threshold`
(the sum of the masses). `length(generator)` is the number of daughters.
The ordinary map reads `generator.random_numbers` coordinates, which is `3n - 4`.

Draw the first event with [`generate_point`](@ref). `rand(generator)` calls
that function. Fixed coordinates use [`generate_from_unit_hypercube`](@ref).

```julia
using RamboOnDiet
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
generator.masses, generator.threshold, generator.random_numbers
```

```julia
using RamboOnDiet, FourVectors
total = FourVector(0.2, 0.0, 0.3; E = 2.1)
generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], total)
generator.total, generator.invariant_mass
```

## Related

- [`generate_point`](@ref) — draw one event from this parent
- [`generate_from_unit_hypercube`](@ref) — map a chosen unit-cube vector
- [`generate_aligned_point`](@ref) — draw an event with the overall rotation fixed
- [`generate_aligned_from_unit_hypercube`](@ref) — map the shorter aligned coordinates
"""
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
    return PhaseSpaceGenerator{T,N,typeof(massesT)}(
        massesT,
        total,
        invariant_mass,
        threshold,
    )
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
