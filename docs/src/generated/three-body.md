```@meta
EditURL = "../../three-body.jl"
```

# Three-body generation and Dalitz processing

This tutorial shows a simple three-body workflow:
1. build the total-state four-vector,
2. generate events with `RemboOnDiet.jl`,
3. move the samples into a `DataFrame`,
4. plot Dalitz histograms with `Plots.jl`.

````@example three-body
using Random
using Statistics
using DataFrames
ENV["GKSwstype"] = "100"
using Plots
using RemboOnDiet
using FourVectors
using ThreeBodyDecays

gr()
default(size = (780, 620), legend = false)
````

## Utility helpers

We will use the invariant-mass pairs `(s23, s12) = (\sigma_1, \sigma_3)` as Dalitz coordinates.

````@example three-body
function events_table(points)
    return DataFrame(
        s23 = [invariant_masses(point.momenta).s23 for point in points],
        s31 = [invariant_masses(point.momenta).s31 for point in points],
        s12 = [invariant_masses(point.momenta).s12 for point in points],
        jacobian = [phase_space_weight(point) for point in points],
    )
end

function dalitz_heatmap(df, ms; bins = 60, weights = nothing, title = "")
    xedges = collect(range(lims1(ms)...; length = bins + 1))
    yedges = collect(range(lims3(ms)...; length = bins + 1))
    xcenters = @. (xedges[1:(end-1)] + xedges[2:end]) / 2
    ycenters = @. (yedges[1:(end-1)] + yedges[2:end]) / 2
    counts = zeros(Float64, bins, bins)

    if isnothing(weights)
        for (x, y) in zip(df.s23, df.s12)
            ix = searchsortedlast(xedges, x)
            iy = searchsortedlast(yedges, y)
            if 1 <= ix <= bins && 1 <= iy <= bins
                counts[iy, ix] += 1
            end
        end
    else
        for (x, y, w) in zip(df.s23, df.s12, weights)
            ix = searchsortedlast(xedges, x)
            iy = searchsortedlast(yedges, y)
            if 1 <= ix <= bins && 1 <= iy <= bins
                counts[iy, ix] += w
            end
        end
    end

    z = Matrix{Float64}(undef, bins, bins)
    for iy = 1:bins, ix = 1:bins
        x = xcenters[ix]
        y = ycenters[iy]
        σs = Invariants(ms; σ1 = x, σ3 = y)
        z[iy, ix] = isphysical(σs, ms) ? counts[iy, ix] : NaN
    end

    boundary = border13(ms; Nx = 500)
    border_x = [point.σ1 for point in boundary]
    border_y = [point.σ3 for point in boundary]

    plt = heatmap(
        xcenters,
        ycenters,
        z;
        xlabel = "m²(Kπ)",
        ylabel = "m²(pK)",
        title = title,
        colorbar_title = isnothing(weights) ? "counts" : "weighted counts",
        aspect_ratio = :equal,
        c = :viridis,
        framestyle = :box,
    )
    plot!(plt, border_x, border_y; color = :white, linewidth = 2)
    return plt
end
````

## Massless three-body phase space

In the massless case, the main sampler gives a constant event Jacobian,
so the unweighted Dalitz histogram is already flat.

````@example three-body
rng = MersenneTwister(7)
sqrt_s = 2.0
massless_masses = [0.0, 0.0, 0.0]
massless_generator = PhaseSpaceGenerator(massless_masses, sqrt_s)
massless_points = [rand(rng, massless_generator) for _ = 1:80_000]
massless_df = events_table(massless_points)
massless_ms = ThreeBodyMasses(0.0, 0.0, 0.0; m0 = sqrt_s)

first(massless_df, 5)

allunique(round.(massless_df.jacobian; digits = 12))

massless_plot = dalitz_heatmap(
    massless_df,
    massless_ms;
    bins = 60,
    title = "Massless 3-body Dalitz population",
)

massless_plot
````

## Massive example: `\Lambda_c^+ \to p K^- \pi^+`

For non-zero masses, the main sampler produces weighted events.
The physically flat phase-space density is recovered by filling the Dalitz histogram
with `weights = jacobian`.

````@example three-body
lc_masses = (m0 = 2.28646, p = 0.93827208816, K = 0.493677, π = 0.13957039)

massive_masses = [lc_masses.p, lc_masses.K, lc_masses.π]
massive_generator = PhaseSpaceGenerator(massive_masses, lc_masses.m0)
massive_points = [rand(rng, massive_generator) for _ = 1:80_000]
massive_df = events_table(massive_points)
massive_ms = ThreeBodyMasses(lc_masses.p, lc_masses.K, lc_masses.π; m0 = lc_masses.m0)

first(massive_df, 5)

std(massive_df.jacobian) / mean(massive_df.jacobian)

massive_plot = dalitz_heatmap(
    massive_df,
    massive_ms;
    bins = 60,
    weights = massive_df.jacobian,
    title = "Lc -> pKπ weighted Dalitz from the main sampler",
)

massive_plot
````

## Building the total four-momentum explicitly

The high-level `generate_point(rng, masses, sqrt_s)` helper works in the center-of-mass frame.
If you want a moving parent, construct the total-state `FourVector` and use `generate_momenta`.

````@example three-body
beam_total = FourVector(0.35, -0.2, 1.1; E = sqrt(lc_masses.m0^2 + 0.35^2 + 0.2^2 + 1.1^2))
boosted_generator = PhaseSpaceGenerator(massive_masses, beam_total)
boosted_point = rand(rng, boosted_generator)

total_momentum(boosted_point.momenta)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

