import ThreeBodyDecays

@testset "ThreeBodyDecays coordinate-map cross-check" begin
    configurations = (
        ((0.0, 0.0, 0.0), 2.0),
        ((0.3, 0.4, 0.5), 2.0),
        ((0.93827208816, 0.493677, 0.13957039), 2.28646),
        ((0.0, 0.2, 0.3), 1.5),
        ((0.3, 0.4, 0.5), 1.201),
    )
    for (masses, M) in configurations
        ms = ThreeBodyDecays.ThreeBodyMasses(masses...; m0 = M)
        generator = PhaseSpaceGenerator(masses, M)
        density = ThreeBodyDecays.phase_space_integrand(_ -> 1.0, ms; k = 1)
        lo, hi = (masses[2] + masses[3])^2, (M - masses[1])^2
        delta = M - sum(masses)
        # Interior mass coordinates avoid singular zero-mass clusters. Polar
        # endpoints are regular and also test the two Dalitz boundary branches.
        for u1 in (0.01, 0.2, 0.6, 0.99), u2 in (0.0, 0.15, 0.5, 0.85, 1.0)
            # Analytic conversion of the INPUTS, independent of generated momenta
            # and the production root solver. x2σs is linear in s23, not in u1.
            v = u1 / (1 + sqrt(1 - u1)) # 1 - sqrt(1 - u1), without cancellation
            m23 = masses[2] + masses[3] + delta * sqrt(v)
            x = ((m23^2 - lo) / (hi - lo), u2)
            reference = ThreeBodyDecays.x2σs(x, ms; k = 1)
            point = generate_aligned_from_unit_hypercube([u1, u2], generator)
            actual = (
                mass2(point[2] + point[3]),
                mass2(point[3] + point[1]),
                mass2(point[1] + point[2]),
            )
            for i = 1:3
                @test actual[i] ≈ reference[i] atol = 2e-11 * M^2 rtol = 2e-10
            end

            # Convert the reference's dx1 dx2 density to du1 du2 and restore
            # our dΦ = δ⁴(P - Σp) ∏d³p/(2E) convention (π²/4 factor).
            dx1_du1 = m23 * delta / (2sqrt(v) * sqrt(1 - u1) * (hi - lo))
            expected_weight = (π^2 / 4) * density(x) * dx1_du1
            @test phase_space_weight(point) ≈ expected_weight rtol = 2e-9
        end
    end
end
