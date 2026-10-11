using LorentzVectorBase: E, px, py, pz, mass, mass2, spatial_magnitude

@testset "Aligned phase space" begin
    rng = MersenneTwister(81)
    @test required_aligned_random_numbers(3) == 2
    @test required_aligned_random_numbers(4) == 5
    @test required_aligned_random_numbers(6) == 11
    @test_throws ArgumentError required_aligned_random_numbers(2)
    @test_throws ArgumentError generate_aligned_point([0.1, 0.2], 1.0)
    generator = PhaseSpaceGenerator([0.2, 0.3, 0.4], 2.0)
    @test_throws ArgumentError generate_aligned_from_unit_hypercube([0.5], generator)
    for bad in (-0.1, 1.1, NaN, Inf)
        @test_throws ArgumentError generate_aligned_from_unit_hypercube(
            [0.5, bad],
            generator,
        )
    end

    function validate_aligned_point(point, masses, total)
        for (p, m) in zip(point, masses)
            @test mass2(p) ≈ m^2 atol = 5e-10
        end
        summed = total_momentum(point.momenta)
        for component in (E, px, py, pz)
            @test component(summed) ≈ component(total) atol = 5e-10
        end
    end

    for n = 3:6, massive in (false, true)
        masses = massive ? collect(range(0.1, 0.4; length = n)) : zeros(n)
        M = sum(masses) + 2.0
        generator = PhaseSpaceGenerator(masses, M)
        moving = FourVector(0.4, -0.3, 0.8; E = sqrt(M^2 + 0.89))
        for _ = 1:25
            rs = rand(rng, required_aligned_random_numbers(n))
            point = generate_aligned_from_unit_hypercube(rs, generator)
            validate_aligned_point(point, masses, generator.total)
            @test px(point[1]) ≈ 0 atol = 1e-12
            @test py(point[1]) ≈ 0 atol = 1e-12
            @test pz(point[1]) < 0
            @test py(point[2]) ≈ 0 atol = 1e-12
            @test px(point[2]) >= 0
            @test isfinite(point.weight) && point.weight > 0
            if !massive
                @test point.weight ≈ massless_phase_space_volume(M^2, n)
            end

            # Check the physical meaning of the first mass coordinate independently.
            cluster = total_momentum(point.momenta[2:end])
            v = ((mass(cluster) - sum(masses[2:end])) / (M - sum(masses)))^2
            @test (n-1) * v^(n-2) - (n-2) * v^(n-1) ≈ rs[1] atol = 1e-10

            # A moving parent must be a single boost of the aligned rest event.
            lab = generate_aligned_from_unit_hypercube(rs, masses, moving)
            validate_aligned_point(lab, masses, moving)
            inverse = FourVector(-px(moving), -py(moving), -pz(moving); E = E(moving))
            for i = 1:n
                recovered = RamboOnDiet.boost(lab[i], inverse)
                for component in (E, px, py, pz)
                    @test component(recovered) ≈ component(point[i]) atol = 5e-10
                end
            end
            @test lab.weight ≈ point.weight
        end
    end

    # Sampler reproducibility and convenience overloads.
    a = generate_aligned_point(MersenneTwister(9), generator)
    b = generate_aligned_from_unit_hypercube(
        rand(MersenneTwister(9), required_aligned_random_numbers(length(generator))),
        generator,
    )
    @test a.weight == b.weight
    @test all(px(a[i]) == px(b[i]) for i = 1:length(a))
    @test length(generate_aligned_point(generator)) == length(generator)
    @test length(generate_aligned_point((0.1, 0.2, 0.3), 2.0)) == 3
    @test length(generate_aligned_point(rng, [0.1, 0.2, 0.3], 2.0)) == 3

    threshold = generate_aligned_from_unit_hypercube([0.4, 0.6], [0.25, 0.25, 0.5], 1.0)
    @test threshold.weight == 0
    @test all(iszero(spatial_magnitude(p)) for p in threshold)

    # Analytic massless 3-body Dalitz moments, including the integrated orientation.
    points = [generate_aligned_point(rng, zeros(3), 2.0) for _ = 1:10_000]
    for pair in ((1, 2), (2, 3), (3, 1))
        xs = [mandelstam(p.momenta, pair...) / 4 for p in points]
        @test mean(xs) ≈ 1/3 atol = 0.01
        @test mean(abs2, xs) ≈ 1/6 atol = 0.01
    end
end
