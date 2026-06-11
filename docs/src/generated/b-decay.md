```@meta
EditURL = "../../b-decay.jl"
```

# Cascade decay: `B^+ \to D^{*+} D^- K^+` with `D^{*+} \to D^0 \pi^+`

This tutorial generates **flat phase-space** events for a sequential decay chain,
stores the four-momenta of the four final-state particles in a `DataFrame`,
and writes `100_000` events to an [Apache Arrow](https://arrow.apache.org/) file.

Decay topology:

```text
B^+  ->  D^{*+}  +  D^-  +  K^+
         |
         +-> D^0 + \pi^+
```

## Importance sampling

`RemboOnDiet.jl` does not return uniformly distributed momenta when masses are non-zero.
Each event comes with a **phase-space weight** (Jacobian factor) such that weighted
histograms and weighted averages reproduce the flat Lorentz-invariant measure.

For this cascade, the full weight factorizes over the sequential steps:

```math
w = w_{B \to D^{*+} D^- K^+} \times w_{D^{*+} \to D^0 \pi^+}.
```

Always use the column `weight` when integrating observables or filling histograms.
Unweighted event samples are biased; the weights correct that bias.

````@example b-decay
using Random
using Statistics
using DataFrames
using Arrow
using RemboOnDiet
using FourVectors
using LorentzVectorBase
````

## Particle masses (GeV)

Masses follow the PDG-style naming used in amplitude analyses.

````@example b-decay
const masses = (
    Bp = 5.27934,
    Dp = 1.86965,       # D^±
    Kp = 0.493677,
    D0 = 1.86483,
    pip = 0.13957039,
    Dst = 2.01026,      # D^{*+}
)
````

The three-body final state is ordered as `(D^{*+}, D^-, K^+)`.

````@example b-decay
const three_body_masses = [masses.Dst, masses.Dp, masses.Kp]
const dst_daughters = [masses.D0, masses.pip]

@assert sum(three_body_masses) < masses.Bp "B mass is below the three-body threshold."
@assert sum(dst_daughters) < masses.Dst "D* mass is below the D0 pi threshold."
````

## Helpers

````@example b-decay
function four_momentum_from_row(row, prefix::AbstractString)
    return FourVector(
        row[Symbol(prefix * "_px")],
        row[Symbol(prefix * "_py")],
        row[Symbol(prefix * "_pz")];
        E = row[Symbol(prefix * "_E")],
    )
end

function momentum_row(prefix::AbstractString, p::FourVector)
    return (
        Symbol(prefix * "_px") => LorentzVectorBase.px(p),
        Symbol(prefix * "_py") => LorentzVectorBase.py(p),
        Symbol(prefix * "_pz") => LorentzVectorBase.pz(p),
        Symbol(prefix * "_E") => LorentzVectorBase.E(p),
    )
end

function cascade_event(rng::AbstractRNG)
    three_body = rand(rng, PhaseSpaceGenerator(three_body_masses, masses.Bp))
    p_Dst, p_Dm, p_Kp = three_body.momenta
    weight_3body = phase_space_weight(three_body)

    dst_decay = generate_momenta(rng, dst_daughters, p_Dst)
    p_D0, p_pip = dst_decay.momenta
    weight_dst = phase_space_weight(dst_decay)

    row = NamedTuple()
    for (name, p) in (
        (:D0, p_D0),
        (:pip, p_pip),
        (:Dm, p_Dm),
        (:Kp, p_Kp),
    )
        row = merge(row, momentum_row(string(name), p))
    end

    return merge(
        row,
        (
            weight_3body = weight_3body,
            weight_dst = weight_dst,
            weight = weight_3body * weight_dst,
        ),
    )
end

function generate_events(rng::AbstractRNG, n_events::Integer)
    rows = Vector{NamedTuple}(undef, n_events)
    for i = 1:n_events
        rows[i] = cascade_event(rng)
    end
    return DataFrame(rows)
end

function validate_kinematics(df::DataFrame; atol = 1e-9)
    p_B = FourVector(0.0, 0.0, 0.0; E = masses.Bp)
    for row in eachrow(df)
        p_D0 = four_momentum_from_row(row, "D0")
        p_pip = four_momentum_from_row(row, "pip")
        p_Dm = four_momentum_from_row(row, "Dm")
        p_Kp = four_momentum_from_row(row, "Kp")

        p_Dst = total_momentum([p_D0, p_pip])
        p_total = total_momentum([p_D0, p_pip, p_Dm, p_Kp])

        @assert isapprox(LorentzVectorBase.mass(p_Dst), masses.Dst; atol, rtol = atol)
        @assert isapprox(LorentzVectorBase.mass(p_total), masses.Bp; atol, rtol = atol)
        @assert isapprox(total_momentum([p_Dst, p_Dm, p_Kp]), p_B; atol, rtol = atol)
        @assert isapprox(p_total, p_B; atol, rtol = atol)
    end
    return nothing
end

dalitz_limits_s23(m, M) = ((m[2] + m[3])^2, (M - m[1])^2)
dalitz_limits_s12(m, M) = ((m[1] + m[2])^2, (M - m[3])^2)

function dalitz_kibble(invs, m, M)
    msq = (m .^ 2)..., M^2
    return kallen(
        kallen(msq[4], msq[1], invs.s23),
        kallen(msq[4], msq[2], invs.s31),
        kallen(msq[4], msq[3], invs.s12),
    )
end

function isphysical_dalitz_point(s23, s12, m, M)
    s31 = M^2 + sum(abs2, m) - s23 - s12
    invs = (s23 = s23, s31 = s31, s12 = s12)
    lo23, hi23 = dalitz_limits_s23(m, M)
    lo12, hi12 = dalitz_limits_s12(m, M)
    return lo23 <= s23 <= hi23 && lo12 <= s12 <= hi12 && dalitz_kibble(invs, m, M) <= 0
end

function dalitz_bin_fraction(xlo, xhi, ylo, yhi, m, M; subdivisions = 5)
    hits = 0
    total = subdivisions^2
    for ix = 1:subdivisions, iy = 1:subdivisions
        x = xlo + (ix - 0.5) * (xhi - xlo) / subdivisions
        y = ylo + (iy - 0.5) * (yhi - ylo) / subdivisions
        hits += isphysical_dalitz_point(x, y, m, M)
    end
    return hits / total
end

function dalitz_bin_counts(s23, s12, weights = nothing; bins = 24)
    m = three_body_masses
    M = masses.Bp
    xedges = collect(range(dalitz_limits_s23(m, M)...; length = bins + 1))
    yedges = collect(range(dalitz_limits_s12(m, M)...; length = bins + 1))
    counts = zeros(Float64, bins, bins)
    for i in eachindex(s23)
        ix = searchsortedlast(xedges, s23[i])
        iy = searchsortedlast(yedges, s12[i])
        if 1 <= ix <= bins && 1 <= iy <= bins
            counts[iy, ix] += isnothing(weights) ? 1.0 : weights[i]
        end
    end
    return counts, xedges, yedges
end

function interior_dalitz_cv(counts, xedges, yedges; bins = 24)
    fractions = [
        dalitz_bin_fraction(
            xedges[ix],
            xedges[ix+1],
            yedges[iy],
            yedges[iy+1],
            three_body_masses,
            masses.Bp;
            subdivisions = 9,
        ) for ix = 1:bins, iy = 1:bins
    ]
    interior = fractions .> 0.95
    values = counts[interior]
    values = values[values .> 0]
    isempty(values) && return Inf
    return std(values) / mean(values)
end

function validate_importance_sampling(df::DataFrame; bins = 24)
    s23 = Float64[]
    s12 = Float64[]
    cosθ_dst = Float64[]
    weights_3body = df.weight_3body
    weights_full = df.weight

    for row in eachrow(df)
        p_D0 = four_momentum_from_row(row, "D0")
        p_pip = four_momentum_from_row(row, "pip")
        p_Dm = four_momentum_from_row(row, "Dm")
        p_Kp = four_momentum_from_row(row, "Kp")
        p_Dst = total_momentum([p_D0, p_pip])

        invs = invariant_masses([p_Dst, p_Dm, p_Kp])
        push!(s23, invs.s23)
        push!(s12, invs.s12)

        p_D0_dst = transform_to_cmf(p_D0, p_Dst)
        p_mag = LorentzVectorBase.spatial_magnitude(p_D0_dst)
        push!(cosθ_dst, LorentzVectorBase.pz(p_D0_dst) / p_mag)
    end

    unweighted_counts, xedges, yedges = dalitz_bin_counts(s23, s12; bins = bins)
    weighted_counts, _, _ = dalitz_bin_counts(s23, s12, weights_3body; bins = bins)
    unweighted_cv = interior_dalitz_cv(unweighted_counts, xedges, yedges; bins = bins)
    weighted_cv = interior_dalitz_cv(weighted_counts, xedges, yedges; bins = bins)

    weight_cv_3body = std(weights_3body) / mean(weights_3body)
    weight_cv_dst = std(df.weight_dst) / mean(df.weight_dst)
    weight_variance_factor = mean(abs2, weights_3body) / mean(weights_3body)^2
    weighted_mean = mean(weighted_counts[weighted_counts .> 0])
    sampling_floor = sqrt(weight_variance_factor * mean(weights_3body) / weighted_mean)
    cosθ_mean = mean(cosθ_dst)
    cosθ_second_moment = mean(abs2, cosθ_dst)

    @assert weight_cv_3body > 0.05 "Expected non-constant three-body importance weights."
    @assert weight_cv_dst < 1e-10 "Two-body D* weight should be constant for fixed masses."
    @assert all(df.weight .≈ df.weight_3body .* df.weight_dst)
    @assert weighted_cv < unweighted_cv / 1.5
    @assert weighted_cv < 3.5 * sampling_floor
    @assert abs(cosθ_mean) < 0.03
    @assert abs(cosθ_second_moment - 1 / 3) < 0.03

    return (
        weight_cv_3body = weight_cv_3body,
        weight_cv_dst = weight_cv_dst,
        dalitz_cv_unweighted = unweighted_cv,
        dalitz_cv_weighted = weighted_cv,
        sampling_floor = sampling_floor,
        cosθ_mean = cosθ_mean,
        mean_weight = mean(weights_full),
    )
end
````

## Generate and inspect a small sample

````@example b-decay
rng = MersenneTwister(2025)
preview = generate_events(rng, 5)
first(preview, 3)
````

## Generate `100_000` phase-space events

````@example b-decay
const N_EVENTS = 100_000
const OUT_DIR = joinpath(@__DIR__, "generated")
const ARROW_FILE = joinpath(OUT_DIR, "b-decay-events.arrow")

mkpath(OUT_DIR)
events = generate_events(rng, N_EVENTS)
validate_kinematics(events)
sampling_report = validate_importance_sampling(events)
sampling_report
````

The checks above confirm:
- the three-body Jacobian varies event-by-event (importance sampling is active),
- the `D^{*+}` two-body Jacobian is constant for fixed masses,
- the weighted three-body Dalitz projection is much flatter than the unweighted one,
- the `D^0` polar angle in the `D^{*+}` rest frame stays isotropic.

````@example b-decay
Arrow.write(
    ARROW_FILE,
    events;
    metadata = Dict(
        "description" => "Flat phase-space events for B+ -> D0 pi+ D- K+ via D*+ cascade",
        "importance_sampling" => "true",
        "weight_column" => "weight",
        "n_events" => string(N_EVENTS),
    ),
)

nrow(events)
ARROW_FILE
isfile(ARROW_FILE)
````

## Read the file back and compute a weighted average

For any observable `f(event)`, the phase-space average is

```math
\langle f \rangle = \frac{\sum_i w_i f_i}{\sum_i w_i}.
```

````@example b-decay
reloaded = DataFrame(Arrow.Table(ARROW_FILE))
nrow(reloaded)

function m2_Dm_Kp(row)
    p_Dm = four_momentum_from_row(row, "Dm")
    p_Kp = four_momentum_from_row(row, "Kp")
    return LorentzVectorBase.mass2(total_momentum([p_Dm, p_Kp]))
end

observable = m2_Dm_Kp.(eachrow(reloaded))
weighted_mean = sum(reloaded.weight .* observable) / sum(reloaded.weight)
unweighted_mean = mean(observable)

weighted_mean
unweighted_mean
````

`weighted_mean` is the phase-space average. `unweighted_mean` is biased and should not
be used for physics integrals.

````@example b-decay
first(reloaded, 3)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

