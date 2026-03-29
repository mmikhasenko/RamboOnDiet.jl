using Statistics
using FHist

function validate_point(point, masses, total; atol, rtol)
    @test length(point) == length(masses)
    for (p, m) in zip(point, masses)
        @test LorentzVectorBase.mass(p) ≈ m atol = atol rtol = rtol
    end
    total_out = total_momentum(point.momenta)
    @test LorentzVectorBase.px(total_out) ≈ LorentzVectorBase.px(total) atol = atol rtol = rtol
    @test LorentzVectorBase.py(total_out) ≈ LorentzVectorBase.py(total) atol = atol rtol = rtol
    @test LorentzVectorBase.pz(total_out) ≈ LorentzVectorBase.pz(total) atol = atol rtol = rtol
    @test LorentzVectorBase.E(total_out) ≈ LorentzVectorBase.E(total) atol = atol rtol = rtol
end

function kibble(invs, masses, M)
    msq = (masses .^ 2)..., M^2
    return kallen(
        kallen(msq[4], msq[1], invs.s23),
        kallen(msq[4], msq[2], invs.s31),
        kallen(msq[4], msq[3], invs.s12),
    )
end

function dalitz_limits_s23(masses, M)
    return (masses[2] + masses[3])^2, (M - masses[1])^2
end

function dalitz_limits_s12(masses, M)
    return (masses[1] + masses[2])^2, (M - masses[3])^2
end

function hist2d_counts(xs, ys, xedges, yedges; weights = nothing)
    if isnothing(weights)
        h = Hist2D((xs, ys); binedges = (xedges, yedges))
    else
        h = Hist2D((xs, ys); binedges = (xedges, yedges), weights = weights)
    end
    return bincounts(h)
end

function flatness_score(counts)
    populated = vec(counts[counts .> 0])
    μ = mean(populated)
    σ = std(populated)
    return (; mean = μ, std = σ, cv = σ / μ, n = length(populated))
end

function reduced_chi2(observed, expected)
    mask = expected .>= 5
    ν = count(mask) - 1
    ν > 0 || return Inf
    residual = observed[mask] .- expected[mask]
    χ2 = sum((residual .^ 2) ./ expected[mask])
    return χ2 / ν
end

function isphysical_dalitz_point(s23, s12, masses, M)
    s31 = M^2 + sum(abs2, masses) - s23 - s12
    invs = (s23 = s23, s31 = s31, s12 = s12)
    lo23, hi23 = dalitz_limits_s23(masses, M)
    lo12, hi12 = dalitz_limits_s12(masses, M)
    return lo23 <= s23 <= hi23 &&
           lo12 <= s12 <= hi12 &&
           kibble(invs, masses, M) <= 0
end

function bin_fraction(xlo, xhi, ylo, yhi, masses, M; subdivisions = 5)
    hits = 0
    total = subdivisions^2
    for ix in 1:subdivisions, iy in 1:subdivisions
        x = xlo + (ix - 0.5) * (xhi - xlo) / subdivisions
        y = ylo + (iy - 0.5) * (yhi - ylo) / subdivisions
        hits += isphysical_dalitz_point(x, y, masses, M)
    end
    return hits / total
end

function validate_threebody_flatness(points, masses, M; bins = 24)
    xs = Float64[]
    ys = Float64[]
    ws = Float64[]
    target = M^2 + sum(abs2, masses)

    for point in points
        invs = invariant_masses(point.momenta)
        @test kibble(invs, masses, M) <= 1e-8
        @test invs.s23 + invs.s31 + invs.s12 ≈ target atol = 1e-10 rtol = 1e-10
        push!(xs, invs.s23)
        push!(ys, invs.s12)
        push!(ws, phase_space_weight(point))
    end

    xedges = collect(range(dalitz_limits_s23(masses, M)...; length = bins + 1))
    yedges = collect(range(dalitz_limits_s12(masses, M)...; length = bins + 1))
    weighted = hist2d_counts(xs, ys, xedges, yedges; weights = ws)
    fractions = [
        bin_fraction(xedges[ix], xedges[ix + 1], yedges[iy], yedges[iy + 1], masses, M; subdivisions = 9)
        for ix in 1:bins, iy in 1:bins
    ]
    interior = fractions .> 0.95
    weighted_stats = flatness_score(weighted[interior])
    weight_variance_factor = mean(abs2, ws) / mean(ws)^2
    sampling_floor = sqrt(weight_variance_factor * mean(ws) / weighted_stats.mean)

    @test weighted_stats.cv < 3.5 * sampling_floor
end

function validate_threebody_massless(points, masses, M; bins = 24)
    ws = [phase_space_weight(point) for point in points]
    @test std(ws) / mean(ws) < 1e-12
    validate_threebody_flatness(points, masses, M; bins = bins)
end

function validate_threebody_massive(points, masses, M; bins = 24)
    ws = [phase_space_weight(point) for point in points]
    @test minimum(ws) ≥ 0
    @test std(ws) / mean(ws) > 1e-3
    validate_threebody_flatness(points, masses, M; bins = bins)
end
