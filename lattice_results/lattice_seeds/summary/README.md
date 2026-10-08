# Six states of the zigzag 6×6 topological lattice

One lattice (zigzag AQH 6 × 6, R = 60 rings, J = 80, φ = π/4, d₂ = 0.0125, κ_ex = 1 on the input
ring 0 and the drop ring 55, N = 64 longitudinal modes, pump alone on the edge supermode σ = 29),
six kinds of nonlinear state as the pump power F₀² and the effective detuning Δ_eff of that
supermode are changed. Each is one entry of the archive `../states/states.npz`; the movies here
were rendered straight from those frames with `../animate_state.py`.

| | state | F₀² | Δ_eff | archive name | movie |
|---|---|---|---|---|---|
| 1 | **below threshold** — linear response, no comb | 500 | 1.00 | `linear_P500_D1.00` | `1_below_threshold.mp4` |
| 2 | **static edge soliton** — one pulse per boundary ring, k = 0 along the edge | 3500 | 4.90 | `soliton1_P3500_D4.90` | `2_static_soliton.mp4` |
| 3 | **nested travelling edge soliton** — one envelope localised in φ and along the edge, circulating | 6500 | 6.60 | `nested_soliton_P6500_D6.60` | `3_nested_soliton.mp4` |
| 4 | **chaotic comb** — developed modulational instability | 6500 | 4.50 | `chaos_P6500_D4.50` | `4_chaotic_comb.mp4` |
| 5 | **soliton crystal** — three bound pulses per boundary ring, k = 0, stationary | 3500 | 4.80 | `turing3_P3500_D4.80` | `5_soliton_crystal.mp4` |
| 6 | **Turing roll** — 14 equidistant pulses per boundary ring, k = 0, stationary | 2500 | 0.10 | `turing14_P2500_D0.10` | `6_turing_roll.mp4` |

Every movie shows the 60 rings on the zigzag geometry; the colour along each ring is the power
|ψ_r(φ)|² inside it (black = none, white = most), the input ring is outlined blue, the drop ring
red. The frame step differs between the movies, because the states move on very different time
scales — it is written in each title, and it matters: a travelling state sampled slower than its
orbit shows a fake slow drift (see "Measurement caveat" in `../STATES.md`).

## 1. Below threshold (F₀² = 500, Δ_eff = 1.00) — `1_below_threshold.mp4`, 0.5 lifetimes per frame

The pump is well below the comb threshold of this lattice (between F₀² = 1000 and 1500). The
light fills the boundary rings uniformly in φ (no pulse, every ring a smooth arc), stays on the
edge — brightest on the short way round from the input to the drop ring, the direction the edge
band carries it — and nothing moves: the frames are identical. All longitudinal modes μ ≠ 0 are empty. This
is the regime of the mini-comb feature map (`LATTICE_EXTENSION.md`): with encoding tones added on
the neighbouring edge supermodes the response is exactly periodic and the fine lines nδ of the
drop port are its features. Lattice power 15.6.

## 2. Static edge soliton (F₀² = 3500, Δ_eff = 4.90) — `2_static_soliton.mp4`, 0.5 lifetimes per frame

One sech-like pulse in φ in every boundary ring, all rings in phase (a k = 0 collective state along
the edge), peak/background ≈ 16, stationary: the pulse sits at the same φ in every frame, and the
power drifts by 10⁻¹⁵ over 800 lifetimes. Reached by a blue→red sweep that ignites chaos, chaos
collapsing onto a locked Turing-3 pattern (3 pulses per ring), then pulse surgery (masking two of
the three). It is one member of a multistable family at this point — CW 12.65, 1 pulse 16.12,
2 pulses 19.80, 3 pulses 23.51, i.e. +3.6 per pulse — and a cold start here falls to CW. Its fine
spectrum is a single line per longitudinal tooth: the edge ladder is NOT lit, so this soliton is
not nested.

## 3. Nested travelling edge soliton (F₀² = 6500, Δ_eff = 6.60) — `3_nested_soliton.mp4`, 0.02 lifetimes per frame

A single envelope about three boundary rings wide, with one sharp φ-pulse inside each of those
rings, running round the 20-ring boundary in one direction (the chirality of the edge band) once
every 0.56 lifetimes — 90 frames are about three laps. It passes the drop port once per lap, so
the drop power is a periodic sawtooth, and its spectrum is a rigid equidistant comb across all
longitudinal teeth: this is the nested (2D) topological comb of Mittal et al. The lap time equals
2π/δ with δ = 10.98 the spacing of the edge supermodes (the mini FSR), i.e. the soliton goes round
the edge "superring" at its own FSR. Reached by a cold fine sweep at this pump: the post-collapse
shelf Δ_eff ∈ [5.6, 6.9] is a staircase 9 → 5 → 1 → 0 envelopes and this is its last step.
Verified 1100 lifetimes. Lattice power 23.7.

## 4. Chaotic comb (F₀² = 6500, Δ_eff = 4.50) — `4_chaotic_comb.mp4`, 0.05 lifetimes per frame

The same pump as the nested soliton, 2.1 closer to resonance: developed modulational-instability
chaos, just before the collapse at Δ_eff ≈ 5.2. Light on every boundary ring and in every
longitudinal mode, the pattern in φ changing from frame to frame without any periodicity, the
lattice power (100, the highest of the four) wandering in time, the fine spectra filled. At this
pump chaos owns the whole blue side (Δ_eff from 0.2 to 5.2): there is no Turing window. Two more
chaotic states are in the archive (`chaos_P3500_D3.10`, `chaos_P6500_D1.00`).

## 5. Soliton crystal (F₀² = 3500, Δ_eff = 4.80) — `5_soliton_crystal.mp4`, 0.5 lifetimes per frame

Three pulses in φ in every boundary ring, at the same three positions in every ring (k = 0 along
the edge, like the static soliton), stationary to machine precision: nothing moves between frames.
The archive calls it "Turing-3", but it is not a Turing roll: the pulses are not equidistant —
they sit at φ/2π = 0, 0.31, 0.88, gaps of 112°, 202° and 45° — their spectrum is the broad sech-like
envelope of a soliton (strongest teeth μ = ±7, ±9; 68 % of the comb beyond |μ| = 6, nothing special
at μ = ±3), and the lattice power is additive, 3.6 per pulse (12.65 → 16.12 → 19.80 → 23.51). It is
a bound state of three solitons, a soliton crystal. This is the locked pattern that the chaos of a blue→red sweep collapses onto at
Δ_eff ≈ 4.1 and that survives to 5.2; the static soliton (2) and the two-pulse state were carved
out of it by pulse surgery, and at this point it coexists with CW (a cold start lands on CW, power
12.65, against 23.51 here — multistability). Its fine spectrum is a single line per longitudinal
tooth, so like the static soliton it is a pattern in φ only, not along the edge.

## 6. Turing roll (F₀² = 2500, Δ_eff = 0.10) — `6_turing_roll.mp4`, 0.5 lifetimes per frame

The pattern that modulational instability makes just above its onset: 14 equidistant pulses in φ in
every boundary ring, at the same 14 positions in every ring (k = 0 along the edge), locked and
stationary — power drift 5 × 10⁻⁷ over 20 lifetimes, no rotation of the pattern. Its spectrum is
what distinguishes a roll from solitons: 99 % of the comb sits on the harmonics μ = ±14, ±28 (the
instability is at μ ≈ 14–16 for this pump), nothing at the other teeth, contrast only 2.4. Reached
by a cold start from noise at these parameters, locked after about 250 lifetimes. Rolls of 14–16
pulses also form at F₀² = 2000–3500 for Δ_eff between 0.05 and 0.2; at 2000 and at Δ_eff = 0.2 they
rotate slowly in φ, at 3500 they stay disordered along the edge for hundreds of lifetimes. This is
the lowest pump of the six states, and the one the pump reaches first when it is raised above the
comb threshold (between 1000 and 1500).

## What each means for a feature map

* **Below threshold**: periodic response under periodic drive — the fine lines nδ are read by a
  Fourier transform of the drop port and are noise-free. The regime used so far (Pendulum −249,
  LunarLander +263).
* **Static soliton**: stationary; the drop power is constant in time and only the longitudinal
  comb carries structure. Untested with encoding tones; the single ring's soliton tolerated only
  ε ≈ 0.02 before it broke.
* **Nested soliton**: exactly periodic with the lap time 2π/δ, so the SAME nδ Fourier read-out
  applies, with a comb that the lattice generates itself. Untested with encoding tones.
* **Soliton crystal**: stationary like the static soliton, three times its contrast in the
  longitudinal comb. Untested with tones.
* **Turing roll**: stationary; the single ring's `rolls` regime needed two-sided tones to pin the
  pattern and tolerated only weak ones (ε = 0.1). Untested here.
* **Chaotic comb**: no periodicity; only a time average of the line powers can be read, and the
  inputs must survive the chaos noise. On the AQH 4 × 4 lattice (F₀² = 800) this regime did not
  learn the pendulum (−889, the level of a linear policy); these 6 × 6 states are pumped harder.

Units: κ_in = κ/2 = 1, time in lifetimes 2/κ, Δ > 0 red-detuned (the opposite sign to the
Topological Photonics Explorer). Rendering: from the repository root,
`python lattice_results/lattice_seeds/animate_state.py --init lattice_results/lattice_seeds/states/states.npz:<archive name> --F2 … --de … --step_lt … --nfr … --out lattice_results/lattice_seeds/summary/<name>`.
