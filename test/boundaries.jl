# A controlled RNG proves that only mass zeros are redrawn, including repeated zeros.
mutable struct BoundaryRNG <: AbstractRNG
    values::Vector{Float64}
    position::Int
end
function Random.rand(rng::BoundaryRNG, ::Type{T}) where {T<:AbstractFloat}
    value = T(rng.values[rng.position])
    rng.position += 1
    return value
end

@testset "Coordinate boundaries" begin
    for n = 3:5, massive in (false, true)
        masses = massive ? fill(0.2, n) : zeros(n)
        generator = PhaseSpaceGenerator(masses, sum(masses) + 2.0)
        for (mapper, count) in (
            (generate_from_unit_hypercube, 3n - 4),
            (generate_aligned_from_unit_hypercube, 3n - 7),
        )
            for i = 1:(n-2), endpoint in (0.0, 1.0)
                rs = fill(0.4, count)
                rs[i] = endpoint
                @test_throws DomainError mapper(rs, generator)
                @test rs[i] == endpoint # deterministic inputs are not clamped
            end
            for endpoint in (0.0, 1.0)
                rs = vcat(fill(0.4, n - 2), fill(endpoint, count - n + 2))
                point = mapper(rs, generator)
                @test all(p -> all(isfinite, p), point)
                @test isfinite(point.weight)
                @test all(isapprox.(sum(point.momenta), generator.total; atol = 1e-10))
            end
            for bad in (-0.1, 1.1, NaN, Inf)
                rs = fill(0.4, count)
                rs[end] = bad
                @test_throws ArgumentError mapper(rs, generator)
            end
        end
    end

    # No mass coordinates exist for two-body events: zero angles are valid.
    @test isfinite(
        generate_from_unit_hypercube([0.0, 1.0], PhaseSpaceGenerator([0.3, 0.4], 2.0)).weight,
    )

    # Preserve the existing global-threshold convention even at cube endpoints.
    threshold = PhaseSpaceGenerator([0.25, 0.25, 0.5], 1.0)
    for endpoint in (0.0, 1.0)
        for (mapper, count) in
            ((generate_from_unit_hypercube, 5), (generate_aligned_from_unit_hypercube, 2))
            point = mapper(fill(endpoint, count), threshold)
            @test point.weight == 0
            @test all(p -> all(isfinite, p), point)
            @test sum(point.momenta) ≈ threshold.total
        end
    end

    for n = 3:5
        generator = PhaseSpaceGenerator(fill(0.2, n), n * 0.2 + 2.0)
        for (sampler, mapper, count) in (
            (generate_point, generate_from_unit_hypercube, 3n - 4),
            (generate_aligned_point, generate_aligned_from_unit_hypercube, 3n - 7),
        )
            draws = vcat(repeat([0.0, 0.0, 0.4], n - 2), zeros(count - n + 2))
            rng = BoundaryRNG(draws, 1)
            actual = sampler(rng, generator)
            expected = mapper(vcat(fill(0.4, n - 2), zeros(count - n + 2)), generator)
            @test actual.momenta == expected.momenta
            @test actual.weight == expected.weight
            @test rng.position == length(draws) + 1
        end
    end
end
