# # Four-body sampling in a three-pion mass window
#
# A 190 GeV/c π− beam on a proton at rest produces `[recoil, π−, π+, π−]`.
# In the four-body map the first unit-cube coordinate is the three-pion cluster
# mass. Drawing that coordinate from the image of `[1.02, 1.7] GeV`, and
# multiplying the weight by the width of the interval, keeps every event inside
# the window.
#
# The weighted sample estimates the restricted measure
#
# ```math
# d\Phi_4
# = \delta^{4}\!\left(P - \sum_i p_i\right)
# \prod_i \frac{d^{3}p_i}{2 E_i},
# ```
#
# without factors of `2π`, the same convention as `phase_space_weight`.
# The events stay weighted. Absolute massive weights in this normalization
# require RamboOnDiet v0.5.0 or later.
#
# The numbers on this page were produced locally from
# `docs/four-body-mass-window.jl` and committed with the manual. Building the
# documentation does not rerun the sample.

# ## Beam and window
#
# The total energy is fixed by the beam momentum. For four bodies the first cube
# input is the three-pion mass through
#
# ```math
# 3v^2 - 2v^3 = u_1,
# \qquad
# m_{3\pi} = 3m_\pi + \Delta\sqrt{v},
# \qquad
# \Delta = W - m_p - 3m_\pi.
# ```
#
# Inverting that map gives the coordinate interval the sampler is allowed to draw.

using RamboOnDiet
using FourVectors
using Random
using LinearAlgebra
using Printf
using Test #hide
struct Page #hide
    body::String #hide
end #hide
function Base.show(io::IO, ::MIME"text/markdown", page::Page) #hide
    print(io, page.body) #hide
end #hide
nothing #hide

const mpi = 0.13957039
const mp = 0.93827208816
const beam_momentum = 190.0
const W = sqrt(mp^2 + mpi^2 + 2mp * sqrt(beam_momentum^2 + mpi^2))
const mass_min, mass_max = 1.02, 1.7
const kinetic_mass = W - mp - 3mpi

function mass_coordinate(m)
    z = ((m - 3mpi) / kinetic_mass)^2
    return z^2 * (3 - 2z)
end

const coordinate_min = mass_coordinate(mass_min)
const coordinate_max = mass_coordinate(mass_max)
const coordinate_width = coordinate_max - coordinate_min
const generator = PhaseSpaceGenerator([mp, mpi, mpi, mpi], W)
note = "The collision energy is **$(@sprintf("%.6f", W)) GeV**. " #hide
note *= "The mass window maps onto `u₁ ∈ [$(@sprintf("%.3e", coordinate_min)), $(@sprintf("%.3e", coordinate_max))]`. " #hide
note *= "That interval is **$(@sprintf("%.6f", 100 * coordinate_width))%** of the cube, " #hide
note *= "so rejection sampling would keep one event in **$(@sprintf("%.0f", 1 / coordinate_width))**." #hide
Page(note) #hide

# ## One event inside the window
#
# Replace the first coordinate by a uniform draw on `[u_a, u_b]` and multiply
# the generator weight by `u_b - u_a`. The second coordinate is the remaining
# mass input and has to stay in the open unit interval. The daughters come out
# as `[recoil, π−, π+, π−]`.

function restricted_event(rng)
    coordinates = rand(rng, 8)
    coordinates[1] = coordinate_min + coordinate_width * coordinates[1]
    ## Mass coordinates are defined on the open unit interval.
    while iszero(coordinates[2])
        coordinates[2] = rand(rng)
    end
    event = generate_from_unit_hypercube(coordinates, generator)
    weight = coordinate_width * phase_space_weight(event)
    return event, weight, coordinates[1]
end
nothing #hide

# ## Reference spectrum
#
# The same restricted volume has an independent factorization. At three-pion
# mass `m` the four-body measure splits into the recoil and the three-pion decay,
#
# ```math
# \frac{d\Phi_4}{dm} = 2m \, \Phi_2(W; m_p, m) \, \Phi_3(m; m_\pi, m_\pi, m_\pi).
# ```
#
# The three-body volume is an integral of two two-body volumes over the pair
# mass squared. The two-body volume below is angle-integrated and uses the
# generator's convention, without factors of `2π`.
# Gauss–Legendre rules of 64 and 128 nodes check that this integral has settled.

phi2(M, a, b) =
    π * sqrt(max(0.0, (M - a - b) * (M + a + b) * (M - a + b) * (M + a - b))) / (2M^2)

function gauss_legendre(n)
    off_diagonal = [k / sqrt(4k^2 - 1) for k = 1:(n-1)]
    decomposition = eigen(SymTridiagonal(zeros(n), off_diagonal))
    return decomposition.values, 2 .* decomposition.vectors[1, :] .^ 2
end

function integrate(f, lower, upper, rule)
    nodes, weights = rule
    midpoint = (lower + upper) / 2
    half_width = (upper - lower) / 2
    return half_width *
           sum(weights[k] * f(midpoint + half_width * nodes[k]) for k in eachindex(nodes))
end

function spectrum(m, rule)
    phi3 = integrate((2mpi)^2, (m - mpi)^2, rule) do t
        phi2(m, mpi, sqrt(t)) * phi2(sqrt(t), mpi, mpi)
    end
    return 2m * phi2(W, mp, m) * phi3
end

function reference_bins(edges, rule)
    return [
        integrate(m -> spectrum(m, rule), edges[b], edges[b+1], rule) for
        b = 1:(length(edges)-1)
    ]
end
nothing #hide

# ## Angle cross-check
#
# Pair masses determine the three-pion decay angles. An independent boost of the
# Cartesian momenta has to reproduce those angles. The comparison uses cyclic
# `(12)3` and `(23)1` conventions and does not enter the weights.

function internal_angles(m3pi, m12, m23)
    s, mu = m3pi^2, mpi^2
    sigma3, sigma1 = m12^2, m23^2
    sigma2 = s + 3mu - sigma1 - sigma3
    kallen(a, b, c) = a^2 + b^2 + c^2 - 2(a * b + a * c + b * c)
    cosine(pair, cross) =
        (2pair * (cross - 2mu) - pair * (s - pair - mu)) /
        sqrt(max(0.0, kallen(pair, mu, mu) * kallen(s, pair, mu)))
    cos_z =
        ((s + mu - sigma1) * (s + mu - sigma3) - 2s * (sigma2 - 2mu)) /
        sqrt(max(0.0, kallen(s, mu, sigma1) * kallen(s, mu, sigma3)))
    return (
        acos(clamp(cosine(sigma3, sigma2), -1, 1)),
        acos(clamp(cosine(sigma1, sigma3), -1, 1)),
        -acos(clamp(cos_z, -1, 1)),
    )
end

spatial3(p) = [p.px, p.py, p.pz]

function rest_spatial(p, frame)
    beta = spatial3(frame) / frame.E
    beta2 = dot(beta, beta)
    beta2 == 0 && return spatial3(p)
    gamma = frame.E / mass(frame)
    return spatial3(p) + ((gamma - 1) * dot(beta, spatial3(p)) / beta2 - gamma * p.E) * beta
end

vector_angle(a, b) = acos(clamp(dot(a, b) / (norm(a) * norm(b)), -1, 1))

function cartesian_angles(p1, p2, p3)
    total = p1 + p2 + p3
    return (
        vector_angle(rest_spatial(p1, p1 + p2), -rest_spatial(p3, p1 + p2)),
        vector_angle(rest_spatial(p2, p2 + p3), -rest_spatial(p1, p2 + p3)),
        -vector_angle(rest_spatial(p1, total), rest_spatial(p3, total)),
    )
end

function angle_discrepancy(momenta, m3pi)
    m12 = mass(momenta[2] + momenta[3])
    m23 = mass(momenta[3] + momenta[4])
    from_masses = internal_angles(m3pi, m12, m23)
    from_boosts = cartesian_angles(momenta[2], momenta[3], momenta[4])
    return maximum(abs.(from_masses .- from_boosts)), all(isfinite, from_masses)
end
nothing #hide

# ## Sample
#
# One hundred thousand events, seed `20261021`. Raw counts follow the restricted
# cube map. The weights estimate the absolute integral of `d\Phi_4` in each mass bin.

function sample_window(; n = 100_000, seed = 20261021)
    rng = MersenneTwister(seed)
    edges = collect(range(mass_min, mass_max; length = 11))
    counts, sums, sums2 = zeros(Int, 10), zeros(10), zeros(10)
    minmass, maxmass = Inf, -Inf
    shell_error = conservation_error = map_error = angle_error = 0.0
    positive_weights = finite_angles = true
    for _ = 1:n
        event, weight, u1 = restricted_event(rng)
        momenta = event.momenta
        parent = sum(momenta)
        m3pi = mass(momenta[2] + momenta[3] + momenta[4])
        discrepancy, finite = angle_discrepancy(momenta, m3pi)
        angle_error = max(angle_error, discrepancy)
        finite_angles &= finite
        positive_weights &= isfinite(weight) && weight > 0
        minmass, maxmass = min(minmass, m3pi), max(maxmass, m3pi)
        map_error = max(map_error, abs(mass_coordinate(m3pi) - u1))
        conservation_error = max(
            conservation_error,
            abs(parent.E - W),
            abs(parent.px),
            abs(parent.py),
            abs(parent.pz),
        )
        shell_error = max(
            shell_error,
            maximum(abs(mass(p) - mi) for (p, mi) in zip(momenta, (mp, mpi, mpi, mpi))),
        )
        bin = clamp(searchsortedlast(edges, m3pi), 1, 10)
        counts[bin] += 1
        sums[bin] += weight
        sums2[bin] += weight^2
    end
    return (;
        n,
        seed,
        edges,
        counts,
        sums,
        sums2,
        minmass,
        maxmass,
        shell_error,
        conservation_error,
        map_error,
        angle_error,
        positive_weights,
        finite_angles,
    )
end
nothing #hide

# ## Result
#
# The table compares the weighted Monte Carlo integral in each bin with the
# 128-point quadrature. Pulls use the Monte Carlo standard error. A second rule,
# with 64 points, checks the quadrature itself.

## @hide-start
grouped(x) = replace(string(round(Int, x)), r"(?<=[0-9])(?=(?:[0-9]{3})+$)" => ",")

function result_page(sample, coarse, exact)
    n = sample.n
    estimates = sample.sums / n
    errors = sqrt.((sample.sums2 / n - estimates .^ 2) / (n - 1))
    pulls = (estimates - exact) ./ errors
    total_mc, total_ref = sum(estimates), sum(exact)
    total_se = sqrt((sum(sample.sums2) / n - total_mc^2) / (n - 1))
    total_pull = (total_mc - total_ref) / total_se
    raw_probability = diff(mass_coordinate.(sample.edges)) / coordinate_width
    raw_pulls =
        (sample.counts - n * raw_probability) ./
        sqrt.(n * raw_probability .* (1 .- raw_probability))
    neff = sum(sample.sums)^2 / sum(sample.sums2)
    refinement = maximum(abs.(coarse ./ exact .- 1))

    @testset "Restricted four-body mass window" begin
        @test mass_min <= sample.minmass <= sample.maxmass <= mass_max
        @test sample.positive_weights
        @test sample.finite_angles
        @test sample.angle_error < 2e-8
        @test sample.conservation_error < 1e-9
        @test sample.shell_error < 1e-9
        @test sample.map_error < 1e-11
        @test refinement < 1e-5
        @test maximum(abs.(raw_pulls)) < 5.5
        @test maximum(abs.(pulls)) < 5.5
        @test abs(total_mc - total_ref) < 5.5total_se
        for m in (mass_min, (mass_min + mass_max) / 2, mass_max), r in (0.1, 0.5, 0.9)
            coordinates = fill(r, 8)
            coordinates[1] = mass_coordinate(m)
            momenta = generate_from_unit_hypercube(coordinates, generator).momenta
            @test mass(momenta[2] + momenta[3] + momenta[4]) ≈ m atol = 1e-9
        end
    end

    lines = String[]
    push!(
        lines,
        "Seed `$(sample.seed)`, **$(grouped(n))** events. Every three-pion mass lies in",
    )
    push!(lines, @sprintf("[%.9f, %.9f] GeV.", sample.minmass, sample.maxmass))
    push!(lines, "")
    push!(
        lines,
        "The restricted volume is **$(@sprintf("%.8g", total_mc)) ± $(@sprintf("%.3g", total_se))**.",
    )
    push!(
        lines,
        "Quadrature gives **$(@sprintf("%.8g", total_ref))** (pull **$(@sprintf("%+.2f", total_pull))**).",
    )
    push!(
        lines,
        "The effective sample size is **$(grouped(neff))** events ($(@sprintf("%.2f", 100neff / n))%).",
    )
    push!(
        lines,
        "The 64-point rule matches the 128-point rule to a relative **$(@sprintf("%.2e", refinement))**.",
    )
    push!(lines, "")
    push!(
        lines,
        "Momentum is conserved to **$(@sprintf("%.2e", sample.conservation_error)) GeV**, masses sit on shell to **$(@sprintf("%.2e", sample.shell_error)) GeV**,",
    )
    push!(
        lines,
        "and the cube map inverts to **$(@sprintf("%.2e", sample.map_error))**. Decay angles agree with the Cartesian boosts to **$(@sprintf("%.2e", sample.angle_error)) rad**.",
    )
    push!(
        lines,
        "Unweighted counts follow the cube map; the largest of those pulls is **$(@sprintf("%.2f", maximum(abs.(raw_pulls))))**.",
    )
    push!(lines, "")
    push!(lines, "| m₃π (GeV) | counts | weighted fraction | reference fraction | pull |")
    push!(lines, "|:--|--:|--:|--:|--:|")
    for b = 1:10
        push!(
            lines,
            @sprintf(
                "| %.3f–%.3f | %d | %.6f | %.6f | %+.2f |",
                sample.edges[b],
                sample.edges[b+1],
                sample.counts[b],
                estimates[b] / total_mc,
                exact[b] / total_ref,
                pulls[b],
            ),
        )
    end
    return Page(join(lines, "\n"))
end
nothing
## @hide-end
sample = sample_window()
coarse = reference_bins(sample.edges, gauss_legendre(64))
exact = reference_bins(sample.edges, gauss_legendre(128))
result_page(sample, coarse, exact) #hide
