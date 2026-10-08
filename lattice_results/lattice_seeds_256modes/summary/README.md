# Five states of the zigzag 6×6 topological lattice at 256 longitudinal modes

Same lattice as `../../lattice_seeds/summary/` (zigzag AQH 6 × 6, R = 60 rings, J = 80, φ = π/4,
d₂ = 0.0125, κ_ex = 1 on the input ring 0 and the drop ring 55, pump alone on the edge supermode
σ = 29), now with **N = 256 longitudinal modes per ring instead of 64**. Each state is one entry of
`../states/states.npz`; the movies were rendered straight from those frames with
`../animate_state.py`.

| | state | F₀² | Δ_eff | archive name | movie |
|---|---|---|---|---|---|
| 1 | **Turing roll** — 15 equidistant pulses per boundary ring, k = 0, static | 2500 | 0.10 | `turing15_P2500_D0.10` | `1_turing_roll.mp4` |
| 2 | **static edge soliton** — one pulse per boundary ring, k = 0 — **metastable only** | 3500 | 4.10 | `metastable_k0_soliton_P3500_D4.10` | `2_static_soliton_metastable.mp4` |
| 3 | **nested travelling edge soliton** — one envelope, localised in φ and along the edge, circulating | 6500 | 5.60 | `nested_soliton_P6500_D5.60` | `3_nested_soliton.mp4` |
| 4 | **chaotic comb** — developed modulational instability | 3500 | 3.00 | `chaos_P3500_D3.00` | `4_chaotic_comb.mp4` |
| 5 | **two nested envelopes** — the last step of the cold sweep at this pump | 6500 | 6.20 | `nested_pair_P6500_D6.20` | `5_nested_pair.mp4` |

Each movie shows the 60 rings on the zigzag geometry. The colour along each ring is the power
|ψ_r(φ)|² inside it (black = none, white = most); the input ring is outlined blue, the drop ring
red. The frame step differs between movies because the states move on very different time scales.
It is written in each title, and it matters: a travelling state sampled slower than its orbit
shows a fake slow drift.

## 1. Turing roll (F₀² = 2500, Δ_eff = 0.10) — 0.5 lifetimes per frame

Fifteen equidistant pulses in φ in every boundary ring, at the same fifteen positions in every
ring (k = 0 along the edge). The light stays on the edge; the second row of rings is faint and the
bulk is dark. The pattern does not rotate: its drift is 0.001 roll periods per lifetime. 100 % of
the comb sits on the harmonics μ = ±15, ±30, … with contrast 3.4. There is one weak oscillation,
at ω ≈ 21 ≈ twice the mini FSR, a beat with the edge supermode two ladder rungs away. It does not
grow (drop-ring peak 3.5 % peak-to-peak, lattice power 0.1 %). Reached by a cold start from
noise; locks after about 120 lifetimes. At 64 modes the same point gave 14 pulses. Lattice power
36.1.

## 2. Static edge soliton (F₀² = 3500, Δ_eff = 4.10) — 0.5 lifetimes per frame — METASTABLE

One sharp pulse in φ in every boundary ring, all at the same φ (k = 0 along the edge), peak 13.4.
Nothing moves across the 45 lifetimes of the movie. **But this is not an attractor at 256 modes.**
The frame comes from the 64-mode static soliton, zero-padded to 256 modes and relaxed for 150
lifetimes; held longer it decays to CW within about 500 lifetimes. No stable static soliton was
found in P = 3500–5500, Δ_eff = 3.6–4.6. At lower pump or redder detuning it decays to CW
(lifetime 20–300 lt, longest near Δ_eff ≈ 4.0–4.1). At higher pump (≥ 4500) or bluer detuning
(≤ 3.7) the edge-uniform pulse breaks up into envelopes that run along the edge. The 64-mode
version (`../../lattice_seeds/summary`, state 2) looked stable to 10⁻¹⁵; that stability was a
truncation effect. Lattice power 18.7.

## 3. Nested travelling edge soliton (F₀² = 6500, Δ_eff = 5.60) — 0.02 lifetimes per frame

A single envelope about three boundary rings wide, with one sharp φ-pulse inside, running round
the 20-ring boundary in one direction once every 0.575 lifetimes (34.8 rings per lifetime), so
90 frames are about three laps. The lap time matches 2π/δ = 0.572 with δ = 10.98 the edge
supermode spacing: the soliton circulates the edge "superring" at its own FSR. Its spectrum is a
rigid equidistant comb across all longitudinal teeth, the nested 2D comb. One envelope in every
sample over 650+ lifetimes. Lattice power 22.9.

**How it was reached differs from 64 modes.** A cold sweep does not reach it here: the
post-collapse shelf at this pump ends with two envelopes (state 5) and drops straight to CW. The
state came from the 64-mode nested soliton, zero-padded to 256 modes and held at a bluer
detuning. At the 64-mode point Δ_eff = 6.60 it decays to CW; at Δ_eff = 5.60 it lives.

## 4. Chaotic comb (F₀² = 3500, Δ_eff = 3.00) — 0.05 lifetimes per frame

Developed modulational-instability chaos. Light on every boundary ring as scattered spots that
change from frame to frame without periodicity, the lattice power wandering between 53 and 55,
the fine spectra filled. At this pump chaos covers Δ_eff ≈ 1.4–3.8 in a cold sweep. Lattice
power 54.3.

## 5. Two nested envelopes (F₀² = 6500, Δ_eff = 6.20) — 0.02 lifetimes per frame

Two envelopes running round the edge, each a few rings wide with one φ-pulse inside, with a
rigid comb across all teeth like state 3. This is what the cold sweep at this pump actually
lands on: the shelf Δ_eff ≈ 5.2–6.22 holds two envelopes and falls to CW at 6.23, with no
single-envelope step in between. Lattice power 24.7.

## What changed relative to 64 modes

| | 64 modes | 256 modes |
|---|---|---|
| Turing roll | 14 pulses, static | 15 pulses, static |
| static soliton | stable (drift 10⁻¹⁵) | metastable, decays within ~500 lt |
| soliton crystal (3 pulses) | stable | decays to CW |
| single nested soliton | cold-sweep accessible at Δ_eff 6.60 | only from a seed, only at bluer Δ_eff (5.60) |
| last step of the 6500 cold sweep | 1 envelope | 2 envelopes |
| lap time of the nested soliton | 0.560 lt | 0.575 lt |

A likely mechanism, not yet tested: at 256 modes d₂μ² reaches ≈ 205, wider than the whole
supermode band (≈ ±160), so the comb's tail can phase-match into other supermodes at
|μ| ≈ 40–70 and leak energy there. 64 modes stop at |μ| = 32 and cut that channel off.

Units: κ_in = κ/2 = 1, time in lifetimes 2/κ, Δ > 0 red-detuned (the opposite sign to the
Topological Photonics Explorer). Rendering, from the repository root:
`python lattice_results/lattice_seeds_256modes/animate_state.py --init lattice_results/lattice_seeds_256modes/states/states.npz:<archive name> --F2 … --de … --step_lt … --nfr 90 --settle 0 --out lattice_results/lattice_seeds_256modes/summary/<name>`.
