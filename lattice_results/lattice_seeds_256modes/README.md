# Lattice seeds at 256 longitudinal modes

This folder repeats the state hunt of `../lattice_seeds/` with **N = 256 longitudinal modes per
ring instead of 64**. The method and parameters are unchanged otherwise: zigzag 6×6 AQH, J = 80,
φ = π/4, d₂ = 0.0125, κ_ex = 1, dt = 0.0025, pump on edge supermode σ = 29, cold blue→red
detuning sweeps, stop-and-hold verification.

The reason for redoing it: a soliton at Δ ≈ 5 has width ≈ √(d₂/Δ) ≈ 0.05 rad, while 64 modes
give a φ-grid of 0.1 rad, so the 64-mode states were under-resolved. **They were.** Several
64-mode results do not survive at 256 modes (see "What changed").

## The states

```bash
python lattice_results/lattice_seeds_256modes/states/load.py --list
python lattice_results/lattice_seeds_256modes/states/load.py --check
python lattice_results/lattice_seeds_256modes/states/load.py nested_soliton_P6500_D5.60 --evolve 30
```

| state | F₀² | Δ_eff | what it is | how it was reached |
|---|---|---|---|---|
| `turing15_P2500_D0.10` | 2500 | 0.10 | **static Turing roll**: 15 pulses per boundary ring, k = 0, no drift | cold start |
| `chaos_P3500_D3.00` | 3500 | 3.00 | **chaotic comb**: power wanders, filled fine spectrum | cold sweep |
| `nested_soliton_P6500_D5.60` | 6500 | 5.60 | **single nested travelling soliton**: one envelope, 0.575 lt per lap, stable 650+ lt | seeded (64-mode state upsampled), not cold-sweep accessible |
| `nested_pair_P6500_D6.20` | 6500 | 6.20 | two-envelope nested travelling state | cold sweep |
| `metastable_k0_soliton_P3500_D4.10` | 3500 | 4.10 | k = 0 single soliton per ring — **metastable only**, decays to CW within ~500 lt | seeded |

**A stable static (k = 0) soliton was not found at 256 modes.** Searched P = 3500–5500,
Δ_eff = 3.6–4.6 from a relaxed 256-mode seed. Lower pump or redder detuning: it decays to CW
(lifetimes 20–300 lt, longest near Δ ≈ 4.0–4.1). Higher pump (≥ 4500) or bluer detuning
(≤ 3.7): the edge-uniform state breaks up into localised travelling envelopes.

## What changed relative to 64 modes

| | 64 modes | 256 modes |
|---|---|---|
| Turing roll (P2500, Δ 0.10) | 14 pulses | 15 pulses, still static |
| P3500 cold sweep | locked crystal shelf Δ 4.3–5.2 | **no shelf**: chaos → narrow step at ≈ 4.0 → CW at 4.1 |
| static k = 0 soliton (P3500, Δ 4.9) | stable, drift 1e-15 | **decays to CW**; nowhere stable in the window searched |
| three-soliton crystal (P3500, Δ 4.8) | stable | decays to CW |
| P6500 post-collapse shelf | Δ 5.2–6.9, staircase 9 → 5 → 1 → 0 | Δ 5.2–6.22, minimum count **2** (2 → 0 at 6.23) |
| single nested soliton | cold-sweep accessible at Δ 6.60 | exists only blue-shifted (Δ 5.60) and only from a seed; at 6.60 it decays |
| circulation period | 0.560 lt | 0.575 lt |

A likely mechanism, **not yet tested**: at 256 modes, d₂μ² reaches ≈ 205, which exceeds the
whole supermode band (≈ ±160). The comb's spectral tail can then phase-match into other
supermodes at μ ≈ 40–70, along the parabolas visible in every dispersion map. 64 modes stop at
|μ| = 32, so that leakage channel was cut off, and states that radiate into it looked stable.

## Folder

- `summary/` — a movie and a 2D portrait of each of the five states, with a README.
- `states/` — `states.npz` (5 frames, 60 × 256 complex each), `states.json`, `load.py`,
  `build_states.py` (provenance).
- Tools: same as `../lattice_seeds/`, with N = 256 and two additions: `sweep.py` saves the
  hold-end field of **every** hold (`frames`), and `verify.py --init npz:frames@<Δ>` starts a hold
  from any sweep point. `surgery.py` takes `--init/--de/--F2/--T`.
- `runs/` (raw sweeps and holds, several GB) is not meant to be committed.
