# # Cascade decay: `B^+ \to D^{*+} D^- K^+` with `D^{*+} \to D^0 \pi^+`
#
# A short end-to-end example: generate flat phase-space events for
#
# ```text
# B^+  ->  D^{*+}  +  D^-  +  K^+
#          |
#          +-> D^0 + \pi^+
# ```
#
# and store `100_000` of them in [Apache Arrow](https://arrow.apache.org/) format.
#
# Non-zero masses require **importance-sampling weights**. For this cascade
#
# ```math
# w = w_{B \to D^{*+} D^- K^+} \times w_{D^{*+} \to D^0 \pi^+},
# ```
#
# and the second factor is constant when `D^{*+}` is on shell. Use the `weight`
# column whenever you integrate observables or fill histograms.

using Random
using Statistics
using DataFrames
using Arrow
using RemboOnDiet
using FourVectors
using LorentzVectorBase

# ## Masses (GeV)

const masses = (
    Bp = 5.27934,
    Dp = 1.86965,
    Kp = 0.493677,
    D0 = 1.86483,
    pip = 0.13957039,
    Dst = 2.01026,
)

three_body_masses = [masses.Dst, masses.Dp, masses.Kp]
dst_weight = phase_space_weight(
    generate_momenta(
        MersenneTwister(0),
        [masses.D0, masses.pip],
        FourVector(0, 0, 0; E = masses.Dst),
    ),
)

# ## Generate `100_000` events

rng = MersenneTwister(2025)
generator = PhaseSpaceGenerator(three_body_masses, masses.Bp)

events = DataFrame(
    let
        b_decay = rand(rng, generator)
        p_Dst, p_Dm, p_Kp = b_decay.momenta
        weight_3body = phase_space_weight(b_decay)

        cosθ, ϕ = 2rand(rng) - 1, 2π * rand(rng)
        p_D0, p_pip = decay_two_body(p_Dst, masses.D0, masses.pip, cosθ, ϕ)

        (
            D0_px = LorentzVectorBase.px(p_D0),
            D0_py = LorentzVectorBase.py(p_D0),
            D0_pz = LorentzVectorBase.pz(p_D0),
            D0_E = LorentzVectorBase.E(p_D0),
            pip_px = LorentzVectorBase.px(p_pip),
            pip_py = LorentzVectorBase.py(p_pip),
            pip_pz = LorentzVectorBase.pz(p_pip),
            pip_E = LorentzVectorBase.E(p_pip),
            Dm_px = LorentzVectorBase.px(p_Dm),
            Dm_py = LorentzVectorBase.py(p_Dm),
            Dm_pz = LorentzVectorBase.pz(p_Dm),
            Dm_E = LorentzVectorBase.E(p_Dm),
            Kp_px = LorentzVectorBase.px(p_Kp),
            Kp_py = LorentzVectorBase.py(p_Kp),
            Kp_pz = LorentzVectorBase.pz(p_Kp),
            Kp_E = LorentzVectorBase.E(p_Kp),
            weight = weight_3body * dst_weight,
        )
    end for _ = 1:100_000
)

first(events, 3)

# ## Write Arrow

const ARROW_FILE = joinpath(@__DIR__, "generated", "b-decay-events.arrow")
mkpath(dirname(ARROW_FILE))
Arrow.write(ARROW_FILE, events)

nrow(events)
ARROW_FILE

# ## Read back and form a weighted average
#
# ```math
# \langle f \rangle = \sum_i w_i f_i \,/ \sum_i w_i
# ```

reloaded = DataFrame(Arrow.Table(ARROW_FILE))

m2_Dm_Kp(row) = LorentzVectorBase.mass2(
    total_momentum((
        FourVector(row.Dm_px, row.Dm_py, row.Dm_pz; E = row.Dm_E),
        FourVector(row.Kp_px, row.Kp_py, row.Kp_pz; E = row.Kp_E),
    )),
)

f = m2_Dm_Kp.(eachrow(reloaded))
sum(reloaded.weight .* f) / sum(reloaded.weight)
