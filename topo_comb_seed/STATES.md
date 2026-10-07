# State catalogue — zigzag 6×6 topological lattice

J = 80, φ = π/4, spin −1 (the default of H_zigzag; the eigenvalues do not depend on it, the chirality does: ring 55 is downstream), d₂ = 0.0125, κ_ex = 1, N = 64 longitudinal modes, R = 60 rings,
pump ring 0, drop ring 55, pumped edge supermode σ = 29 (λ = −5.24). Units κ_in = κ/2 = 1,
time in lifetimes 2/κ; Δ = (2/κ)(ω_res − ω_p), so **Δ > 0 is red-detuned** (the opposite sign
convention to the Topological Photonics Explorer, where δ = −Δ).

Three pump powers are kept: **500** (below threshold, linear reference), **3500** (the
stationary k = 0 soliton family) and **6500** (the staircase down to the single nested
travelling soliton). Every state below is in `states/states.npz` under the name in the first
column — load it with `states/load.py`, no raw run data needed.

| archive name | F₀² | Δ_eff | signature | how it was prepared | status |
|---|---|---|---|---|---|
| `linear_P500_D1.00` | 500 | 1.00 | φ-uniform, hold mean = hold end everywhere, dispersion map is the pure noise floor; the tilted resonance peaks at Δ ≈ 0.7 | cold sweep | verified (below threshold; the comb threshold lies between F₀² = 1000 and 1500) |
| `minicomb_P3500_D1.40` | 3500 | 1.40 | discrete equidistant fine lines in tooth μ = 4 under a smooth envelope, inside the dynamic zone | cold sweep | **candidate** — nesting not established; needs stop-and-hold + ladder-vs-rigid-grid check |
| `chaos_P3500_D3.10` | 3500 | 3.10 | hold-end jitter, filled fine spectra, speckled spatiotemporal map | cold sweep | verified (dynamics) |
| `cw_P3500_D4.80` | 3500 | 4.80 | φ-uniform, lattice power 12.65, single fine line | cold start at this detuning | verified |
| `turing3_P3500_D4.80` | 3500 | 4.80 | 3 bound φ-pulses at 0°, 112° and 315° (NOT equidistant: a soliton crystal, not a Turing roll — soliton-type spectrum, power additive 3.6 per pulse), every boundary ring in phase (k = 0); stationary to machine precision; power 23.51 | warm start — a cold start at the same point falls to CW instead (**multistability**) | verified 500 lt |
| `soliton2_P3500_D4.90` | 3500 | 4.90 | 2 pulses per ring, power 19.80 | pulse surgery on turing3 (`surgery.py 2`) | verified 500 lt |
| `soliton1_P3500_D4.90` | 3500 | 4.90 | **1 sech pulse per boundary ring**, edge-rigid (k = 0, **not nested**); power 16.12, peak/background ≈ 16, drift 1e-15 | pulse surgery on turing3 (`surgery.py 1`) | verified 800+ lt; portrait + formation movie |
| `chaos_P6500_D1.00` | 6500 | 1.00 | power wanders, spectra filled — chaos already owns the blue side at this pump | cold sweep | verified 300 lt |
| `chaos_P6500_D4.50` | 6500 | 4.50 | developed chaos, just before the collapse at 5.2 | cold sweep | verified 150 lt |
| `gas_P6500_D5.50` | 6500 | 5.50 | ~10 pulses scattered over the boundary × φ plane; **metastable** — evaporated one envelope at t ≈ 250 lt (power 32 → 30.2) | cold sweep onto the post-collapse shelf | verified 300 lt |
| `gas_P6500_D5.90` | 6500 | 5.90 | ~9 pulses, disordered, localised in both ring and φ; dense fine comb | cold sweep | verified 300 lt |
| `crystal_P6500_D6.30` | 6500 | 6.30 | ~5 pulses, sparse and 2D-localised; clean periodic power sawtooth; rigid equidistant fine comb across all μ | cold sweep | verified 300 lt |
| **`nested_soliton_P6500_D6.60`** | 6500 | 6.60 | **ONE 2D-localised envelope** — ~3 boundary rings wide, a single φ-pulse inside each, contrast 15 — **circulating the 20-ring edge once every 0.560 lifetimes** (35.7 rings/lt, one chirality, zero φ-drift); rigid equidistant comb across all μ; one drop-port pass per lap | cold sweep, no surgery | **verified 1100 lt** |
| `cw_P6500_D6.90` | 6500 | 6.90 | φ-uniform, the shelf is over | cold sweep | verified 300 lt |

Power quantisation at 3500, Δ ≈ 4.9: 12.65 (CW) + 3.6 per pulse → 16.12, 19.80, 23.51. All
four members co-exist at the same parameters; which one you get depends entirely on history.

## The two staircases

**F₀² = 3500** — a cold ramp ignites chaos (Δ ≈ 0.2–3.6), chaos collapses in discrete steps
(3.6–4.3) onto a long locked shelf at power ≈ 27.5 (Δ ∈ [4.3, 5.2]), which dies straight to CW
at 5.2 **with no further sub-steps**. The shelf carries the three-soliton crystal ("turing3" in the archive, a misnomer); the 2- and 1-pulse members
of the family are reachable only by surgery or a warm start, never by the cold sweep itself.

**F₀² = 6500** — chaos from Δ ≈ 0.2 to 5.2 (it owns the blue side too: there is **no Turing
roll at this pump**), then the post-collapse shelf Δ ∈ [5.6, 6.9] resolves under a fine sweep
(step 0.025) into a soliton staircase: **~9 envelopes (5.9) → ~5 (6.3) → 1 (6.6) → 0 (6.9)**.
The last step before CW is the single nested travelling soliton, and the cold sweep lands on
it by itself.

Why it was not found at higher pump: at F₀² = 8000 the same staircase bottoms out at **2**
envelopes (a corner-pinned standing pair), and twelve surgical routes there all failed — hard
carves collapse to CW, carving one of a pair regrows the partner, adiabatic walks in detuning
or in pump die rather than shed an envelope, and a diagonal walk of the 3500 k = 0 soliton
terminates near F₀² ≈ 4700. The minimum envelope count is a property of the pump power; only
at 6500 is it 1.

## Measurement caveat worth remembering

The circulation period is **0.560 lt**, not the ~61 lt a stop-and-hold power plot suggests:
that plot samples every 5 lt, far slower than the orbit, so its apparent slow beat is an
aliasing artifact. Track the envelope centroid with ≥ 10 samples per lap (0.01 lt sampling
gives 56) before quoting any travelling-state period — and use the same step for movie frames.

## Still open

- **Nested / boundary-periodic Turing** — a *regular* k ≠ 0 pattern along the edge, as opposed
  to the disordered gas: not seen at 500/3500/6500. Edge-uniform rolls live at lower pump (an
  8-pulse k = 0 roll was seen at F₀² = 2000, Δ ≈ 3.2, before that data was dropped).
- **Breathers** — no zone found yet where the hold end oscillates regularly rather than
  chaotically at fixed detuning.
- `minicomb_P3500_D1.40` needs the ladder-vs-rigid-grid check to decide whether it is nested.
