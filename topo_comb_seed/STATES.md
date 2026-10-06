# State catalogue — zigzag 6×6 (runs/ is organised per pump power: runs/P3500, P8000, P20000 (P = pump power F0^2)), J = 80, φ = π/4, d₂ = 0.0125, κ_ex = 1 (units κ/2 = 1)

Every distinct state this system has shown us. One row per state: where it lives, how to
recognise it, and the exact seed that re-prepares it. Keep this current — any new state from
any run gets a row, interesting or embarrassing.

| # | state | parameters (F₀², Δ_eff) | signature | seed / preparation | status |
|---|---|---|---|---|---|
| 1 | **CW** (flat lower branch) | any F₀², red of resonance | φ-uniform; lattice P = 12.65 @3500/4.8; single fine line | cold start (noise) anywhere red | verified |
| 2 | **MI chaotic comb** | F₀² ≳ 1500; Δ_eff ≈ 1–4 @3500 | hold-end jitter; filled fine spectra; speckled spatiotemporal | sweep snapshot `P3500/sweep_P3500_J80.npz:snap_+3.10` | verified (dynamics); λ_max not yet measured here |
| 3 | **Locked mini-comb** (equidistant fine lines, smooth envelope; possibly nested Turing — nesting NOT established) | 3500, Δ_eff ≈ +1.4 (also seen at 2000, +1.4, first round) | discrete equidistant fine lines in tooth μ=4 with smooth envelope; state within the "dynamic" zone | `P3500/sweep_P3500_J80.npz:snap_+1.40` | **candidate — needs stop-and-hold + ladder-vs-rigid-grid check** |
| 4 | **Edge-uniform Turing-3** (3 pulses/ring, k = 0 along boundary) | 3500, Δ_eff ∈ [≈4.2, 5.20]; dies 5.25 | stationary (machine precision); 3 equidistant φ-pulses, all edge rings in phase; single fine line | `P3500/verify_P3500_D+4.80_warm.npz:trace` (warm only — cold start falls to #1) | verified 500 lt |
| 5 | **Edge-uniform Turing-8** | 2000, Δ_eff ∈ [3.15, 3.30] | as #4 with 8 pulses, contrast 4.5 | first-round data deleted; re-prepare: sweep at 2000, snap +3.10, warm verify | verified (first round), seed to regenerate |
| 6 | **Two-soliton state** | 3500, 4.9 | 2 pulses; P = 19.80 = CW + 2×3.6 | `P3500/surgery_keep2.npz:a` | verified 500 lt |
| 7 | **Single stationary edge soliton** (HERO; **NOT nested** — uniform along the boundary, k = 0) | 3500, 4.9 | 1 sech pulse/ring, edge-rigid; P = 16.12; drift 1e-15; peak/bg ≈ 16 | `P3500/surgery_keep1.npz:a` (surgery on #4) | verified 800+ lt; portraits + animation |
| 8 | power-quantisation family #1/4/6/7 | 3500, 4.9 | P = 12.65 + 3.6 × n_pulses, n = 0..3, co-existing | — | verified |

| 9 | **Travelling nested comb** (locked envelope train, ~5-6 envelopes circulating the boundary) | 8000, Δ_eff ∈ ≈[6, 8]; verified at 7.0 | lattice power ~flat, DROP power oscillates periodically (circulation ~0.9 lt); staggered blobs along the unrolled boundary (k ≠ 0); rigid equidistant fine comb (spacing ≈ 7) across ALL μ — straight vertical lines leaving the cold guides | `P8000/sweep_P8000_J80.npz:snap_+7.00` → `P8000/verify_P8000_D+7.00_nested.npz` | **verified 300 lt** |
| 10 | second rising branch | 8000, Δ_eff ∈ [8, 12+] | power climbs again with moderate jitter; μ=4 at +9.00 shows few clean lines | `P8000/sweep_P8000_J80.npz:snap_+9.00` | candidate, unexplored |
| 11 | deep chaos sea | 20000, everywhere in [−3, 25] | end-jitter everywhere, filled maps, no locked window | `P20000/sweep_P20000_J80.npz` snaps | verified (negative result) |
| 12 | **Two-envelope travelling state** (the clean minimum of the train family) | 8000, Δ_eff ∈ [7.8, ≈8.05] (50-lt holds; dies 8.06) | 2 envelopes, boundary contrast 23–30; drop oscillates; **keep-1 surgery REGROWS the partner → 2 is the attractor** | `P8000/trainwalk.npz:state_+7.90`, `P8000/microwalk.npz:state_+8.02` | verified |
| 13 | **Edge soliton gas** (~9 pulses, 2D-localised, staggered in ring AND φ) | 6500, Δ_eff ≈ 5.9 (shelf Δ∈[5.6, 6.9] after fine-sweep collapse) | power 29.4 with small regular oscillation; discrete bright spots scattered over the unrolled-boundary × φ plane; dense fine comb | `P6500/sweep_P6500_J80_fine.npz:snap_+5.90` → `P6500/verify_P6500_D+5.90_shelf59.npz` | verified 300 lt; at Δ=5.5 the denser gas (~11) is metastable — evaporated one envelope at t≈250 lt (32→30.2) |
| 14 | **Nested soliton crystal / few-envelope state** (~5 pulses, 2D-localised) | 6500, Δ_eff ≈ 6.3 | power 25.55 with CLEAN periodic sawtooth (~15 lt); sparse 2D-localised pulses; **rigid equidistant fine comb across ALL μ (nested signature)** | `P6500/sweep_P6500_J80_fine.npz:snap_+6.30` → `P6500/verify_P6500_D+6.30_shelf63.npz` | verified 300 lt |
| 15 | **SINGLE NESTED TRAVELLING EDGE SOLITON** (THE TARGET: one envelope, localised in boundary-ring AND φ, circulating the edge) | 6500, Δ_eff = 6.60 (last step before CW at 6.9; shelf staircase 9→5→1→0) | ONE 2D-localised envelope (~3 rings wide, single φ-pulse, contrast 15); circulates the 20-ring boundary in ~61 lt (0.33 rings/lt, one chirality, φ-drift 0); power 23.4 with clean periodic oscillation (one drop-port pass per lap); rigid equidistant fine comb across ALL μ; clean discrete μ=4 comb with smooth envelope | `P6500/sweep_P6500_J80_fine.npz:snap_+6.60` → `P6500/verify_P6500_D+6.60_shelf66.npz:trace` | **verified 300 lt; COLD-SWEEP ACCESSIBLE** |

## Wanted (not yet found)

| state | expected signature | hunt plan |
|---|---|---|
| ~~**SINGLE travelling nested soliton**~~ | **FOUND — state #15** (P6500, Δ 6.60, via the fine-sweep staircase; no surgery needed, cold sweep lands on it) | twelve surgical routes failed at P8000 before the fine 6500 sweep handed it over: | hard carve at 7.0/7.6 (collapse), carve from pair at 7.9 (partner regrows), adiabatic walks (die 8.06, no 2→1), heal+restore (→3), power-walk down (dies 7400), σ=34 end-pumping (no locked branch), intermediate pumps P5000/P6500 (no travelling branch at all — post-collapse is pure CW; the "X" structure in the P6500 +7.50 map is a row-norm transient artifact, verified CW), red-edge carves at 8.05/8.10 (collapse to CW — branch end, not pair interaction), diagonal walk of the P3500 hero (k=0 single-pulse branch terminates ≈P4700, Δ5.6, via a growing breather; cannot reach 8000). In flight: chiral phase-kick carves (select one travelling direction). Next knobs: larger lattice; different pump supermode |
| **Nested / boundary-periodic Turing** (pattern along the edge) | few ladder lines lit and locked; k ≠ 0 structure along unrolled boundary | same sweeps, blue-zone snapshots |
| breathers (if any) | periodic hold-end oscillation at fixed detuning, single frequency in ESA | flag any zone where end-vs-mean oscillates regularly, not chaotically |

Caveat log: #3 was spotted in the first-round summaries (deleted data); the 3500 rerun
regenerates its snapshot. "Nested" claims require the fine-spectrum ladder check — a locked
mini-comb per ring is necessary but not sufficient.

DATA PURGE (2026-10-05 evening, user request): runs/ now keeps ONLY P500, P3500, P6500.
Seeds referenced below under P2000/P8000/P20000 are deleted — regenerate via the quoted
sweep settings if needed (states #5, #9–#12 and the P8000 surgery/walk chain).

P6500 state census (verified): NO regular Turing roll at this pump — chaos already owns the
blue side at Δ=1.0 (verify_P6500_D+1.00), chaos confirmed at 4.5, and the post-collapse shelf
is a disordered soliton gas, not a periodic roll. Edge-uniform Turing rolls live at lower
pump (#5 at P2000, #4 at P3500).

## Waterfall survey (2026-10-05, runs/pumpmode_waterfall.png + combpower_waterfall.png)

Cold blue->red sweeps P = 500..8000 step 500, recording the pumped-supermode mu=0 tooth
|b_{sigma_p,0}|^2 and comb power (total - pump tooth); per-F data in runs/P*/pumpmode_P*.npz.
- Comb threshold refined: between P = 1000 (no comb anywhere) and 1500 (clean locked
  rectangular comb window, Delta_eff 1.5-2.5).
- Chaos collapse front marches red roughly linearly: ~1.5 at P1500 -> ~5.7 at P8000.
- Red edge of the plateau shows discrete state-switching steps (P3000-4000); the P3500
  shelf at 4-5.2 is the Turing-3/soliton family step.
- The travelling-train branch (#9) appears in cold sweeps as a broad comb-power bump at
  Delta 6-8, present at P7500/8000, nearly gone at P7000: cold-accessible birth at
  P ~ 7000-7500 (independent confirmation of the power-walk death at 7400).
- Every panel's far-red rising tail (even below threshold) is LINEAR excitation of the next
  edge-supermode rung (mini-FSR ~10.6), not comb.
