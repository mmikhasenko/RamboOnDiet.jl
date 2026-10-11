# Independent references in the convention dΦ = δ⁴(P - Σp) ∏d³p/(2E).
function reference_two_body(A, b, c)
    q = sqrt((A - b - c) * (A + b + c) * (A - b + c) * (A + b - c)) / (2A)
    return π * q / A
end

# Differentiate the sequential mass map directly, without the production
# massless-volume helper, density ratios, or root solver.
function reference_weight(rs, masses, M)
    n = length(masses)
    K = M - sum(masses)
    clusters = [M]
    jacobian = 1.0
    for i = 2:(n-1)
        a = n - i
        lo, hi = 0.0, 1.0
        for _ = 1:70
            x = (lo + hi) / 2
            if (a + 1) * x^a - a * x^(a + 1) < rs[i-1]
                lo = x
            else
                hi = x
            end
        end
        x = (lo + hi) / 2
        C = K * sqrt(x) + sum(masses[i:end])
        jacobian *= C * K / (sqrt(x) * a * (a + 1) * x^(a - 1) * (1 - x))
        K *= sqrt(x)
        push!(clusters, C)
    end
    push!(clusters, masses[end])
    return jacobian *
           prod(reference_two_body(clusters[i-1], masses[i-1], clusters[i]) for i = 2:n)
end

midpoint_integral(f, a, b, N) =
    (b - a) / N * sum(f(a + (i - 0.5) * (b - a) / N) for i = 1:N)

function integrated_three_body_weight(masses, M, N)
    generator = PhaseSpaceGenerator(masses, M)
    # Smooth the endpoint singularities with r = sin²(πt/2).
    return midpoint_integral(0.0, 1.0, N) do t
        rs = [sinpi(t / 2)^2, 0.3, 0.3, 0.3, 0.3]
        phase_space_weight(generate_from_unit_hypercube(rs, generator)) * π * sinpi(t) / 2
    end
end

@testset "Absolute phase-space normalization" begin
    @testset "Two-body convention" begin
        for masses in ([0.0, 0.0], [0.3, 0.5])
            generator = PhaseSpaceGenerator(masses, 2.0)
            weight = phase_space_weight(generate_from_unit_hypercube([0.3, 0.7], generator))
            @test weight ≈ reference_two_body(2.0, masses...) rtol = 1e-12
            if all(iszero, masses)
                @test weight ≈ π / 2 rtol = 1e-12
            end
        end
    end

    @testset "Sequential-decay Jacobian" begin
        for n = 2:6, r in (0.15, 0.3, 0.75)
            masses = collect(range(0.1, 0.4; length = n))
            M = sum(masses) + 1.3
            rs = collect(range(r, 0.8; length = 3n - 4))
            generator = PhaseSpaceGenerator(masses, M)
            weight = phase_space_weight(generate_from_unit_hypercube(rs, generator))
            @test weight ≈ reference_weight(rs, masses, M) rtol = 1e-10
        end
    end

    @testset "Integrated massive three-body volume" begin
        configurations = (
            ([0.3, 0.4, 0.5], 2.0),
            (fill(0.13957039, 3), 1.86965),
            ([0.93827208816, 0.493677, 0.13957039], 2.28646),
        )
        for (masses, M) in configurations
            reference = midpoint_integral(sum(masses[2:3]), M - masses[1], 20_000) do C
                2C *
                reference_two_body(M, masses[1], C) *
                reference_two_body(C, masses[2], masses[3])
            end
            coarse = integrated_three_body_weight(masses, M, 10_000)
            fine = integrated_three_body_weight(masses, M, 20_000)
            # Midpoint error at square-root endpoints is O(N^(-3/2));
            # refinement also checks that cube quadrature error is below tolerance.
            @test coarse ≈ fine rtol = 2e-6
            @test fine ≈ reference rtol = 2e-6
        end
    end

    @testset "Positive-mass threshold scaling" begin
        for n = 2:6
            masses = fill(0.2, n)
            rs = fill(0.3, 3n - 4)
            weights = map((1e-4, 2e-4)) do delta
                generator = PhaseSpaceGenerator(masses, sum(masses) + delta)
                phase_space_weight(generate_from_unit_hypercube(rs, generator))
            end
            # At fixed interior coordinates the leading Δ power is the same
            # as for the volume; these gaps are above the threshold shortcut.
            @test weights[2] / weights[1] ≈ 2.0^((3n - 5) / 2) rtol = 2e-3
        end
    end
end
