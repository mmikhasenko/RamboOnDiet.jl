using Test
using Random
using RemboOnDiet
using FourVectors
using LorentzVectorBase
using StaticArrays

include("validation.jl")

const JULIA_SEED = 20260329

@testset "RemboOnDiet" begin
    rng = MersenneTwister(JULIA_SEED)

    @testset "Random number bookkeeping" begin
        @test required_random_numbers(2) == 2
        @test required_random_numbers(3) == 5
        @test required_random_numbers(6) == 14
    end

    @testset "Generator object API" begin
        generator = PhaseSpaceGenerator([0.3, 0.4, 0.5], 2.0)
        point = rand(rng, generator)
        validate_point(point, generator.masses, generator.total; atol = 5e-11, rtol = 5e-11)
        @test generator.random_numbers == required_random_numbers(length(generator.masses))
        @test generator.threshold == sum(generator.masses)
        @test generator.invariant_mass ≈ 2.0 atol = 1e-12 rtol = 1e-12

        static_generator = PhaseSpaceGenerator(SVector(0.3, 0.4, 0.5), 2.0)
        @test static_generator.masses isa SVector{3,Float64}
        @test static_generator.random_numbers == 5
    end

    @testset "Massless constant weights" begin
        sqrt_s = 5.0
        for n in 2:6
            if n == 3
                generator = PhaseSpaceGenerator(zeros(n), sqrt_s)
                weights = [phase_space_weight(generate_point(rng, generator)) for _ in 1:64]
                @test std(weights) / mean(weights) < 1e-12
            else
                point = generate_point(rng, zeros(n), sqrt_s)
                expected = float(massless_phase_space_volume(sqrt_s^2, n))
                @test phase_space_weight(point) ≈ expected rtol = 1e-12
            end
        end
    end

    @testset "Conservation and on-shell constraints" begin
        for n in 2:6
            masses = rand(rng, n) .* 0.6
            sqrt_s = sum(masses) + 1.5
            total = FourVector(0.7, -0.4, 0.9; E = sqrt(sqrt_s^2 + 0.7^2 + 0.4^2 + 0.9^2))
            for _ in 1:60
                point = generate_momenta(rng, masses, total)
                validate_point(point, masses, total; atol = 5e-11, rtol = 5e-11)
                @test phase_space_weight(point) ≥ 0
            end
        end
    end

    @testset "Two-body isotropy in the center-of-mass frame" begin
        masses = [0.7, 1.1]
        sqrt_s = 4.2
        cosines = Float64[]
        for _ in 1:3_000
            point = generate_point(rng, masses, sqrt_s)
            push!(
                cosines,
                LorentzVectorBase.pz(point[1]) / LorentzVectorBase.spatial_magnitude(point[1]),
            )
        end
        @test abs(mean(cosines)) < 0.03
        @test abs(mean(abs2, cosines) - 1 / 3) < 0.03
    end

    @testset "Three-body Dalitz checks through the main sampler" begin
        massless_generator = PhaseSpaceGenerator(zeros(3), 2.0)
        massless_events = [generate_point(rng, massless_generator) for _ in 1:12_000]
        validate_threebody_massless(massless_events, massless_generator.masses, 2.0)

        massive_masses = [0.93827208816, 0.493677, 0.13957039]
        massive_generator = PhaseSpaceGenerator(massive_masses, 2.28646)
        massive_events = [generate_point(rng, massive_generator) for _ in 1:12_000]
        validate_threebody_massive(massive_events, massive_masses, 2.28646)
    end

    @testset "Boosted three-body generation" begin
        masses = [0.93827208816, 0.493677, 0.13957039]
        m0 = 2.28646
        total = FourVector(0.4, -0.3, 1.2; E = sqrt(m0^2 + 0.4^2 + 0.3^2 + 1.2^2))
        generator = PhaseSpaceGenerator(masses, total; dalitz_grid = 4097)
        for _ in 1:150
            point = generate_point(rng, generator)
            validate_point(point, masses, total; atol = 5e-11, rtol = 5e-11)
        end
    end
end
