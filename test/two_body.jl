using LorentzVectorBase: E, px, py, pz, mass2, spatial_magnitude

@testset "Explicit-angle two-body decay" begin
    for masses in ([0.0, 0.0], [0.3, 0.5]), moving in (false, true)
        M = 2.0
        total = moving ? FourVector(0.4, -0.3, 0.8; E = sqrt(M^2 + 0.89)) :
                FourVector(0.0, 0.0, 0.0; E = M)
        generator = PhaseSpaceGenerator(masses, total)
        for c in (-1.0, 0.4, 1.0), phi in (0.0, π/3, 2π, -π/3)
            point = decay_two_body(generator; cosθ = c, ϕ = phi)
            expected = generate_from_unit_hypercube([(c + 1)/2, mod(phi, 2π)/2π], generator)
            for i = 1:2
                @test mass2(point[i]) ≈ masses[i]^2 atol = 1e-12
            end
            for component in (E, px, py, pz)
                @test component(sum(point.momenta)) ≈ component(total) atol = 1e-12
            end
            @test point.weight ≈ reference_two_body(M, masses...)
            for component in (E, px, py, pz), i = 1:2
                @test component(point[i]) ≈ component(expected[i]) atol = 1e-12
            end
            if !moving
                q = spatial_magnitude(point[1])
                @test pz(point[1]) ≈ q*c atol = 1e-12
                @test px(point[1]) ≈ q*sqrt(1-c^2)*cos(phi) atol = 1e-12
                @test py(point[1]) ≈ q*sqrt(1-c^2)*sin(phi) atol = 1e-12
            end
        end
    end
    for total in (FourVector(0.0, 0.0, 0.0; E = 1.0),
                  FourVector(0.0, 0.0, 0.75; E = 1.25))
        point = decay_two_body(PhaseSpaceGenerator([0.4, 0.6], total); cosθ = -1, ϕ = 0)
        validate_point(point, [0.4, 0.6], total; atol = 5e-11, rtol = 5e-11)
        @test point.weight == 0
    end
    generator = PhaseSpaceGenerator([0.3, 0.5], 2.0)
    for c in (-1.1, 1.1, NaN, Inf)
        @test_throws ArgumentError decay_two_body(generator; cosθ = c, ϕ = 0)
    end
    for phi in (NaN, Inf, -Inf)
        @test_throws ArgumentError decay_two_body(generator; cosθ = 0, ϕ = phi)
    end
    @test_throws ArgumentError decay_two_body(PhaseSpaceGenerator(zeros(3), 2.0); cosθ = 0, ϕ = 0)
    @test_throws ArgumentError decay_two_body(PhaseSpaceGenerator([0.4, 0.6], 0.9); cosθ = 0, ϕ = 0)
end
